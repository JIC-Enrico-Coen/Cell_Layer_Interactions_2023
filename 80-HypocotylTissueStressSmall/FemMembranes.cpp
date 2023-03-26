#include "FemMembranes.hpp"
#include <MeshBuilder.hpp>
#include <CCDivideCell.hpp>
#include <Geometry.hpp>

namespace mdx
{
  // Convert from Green strain
  double fromGreen(double s)
  {
    if(s == 0)
      return 0;
    else if(s > 0)
      return sqrt(s);
    else
      return -sqrt(-s);
  }

  typedef AttrMap<CCIndex, DoubleVec> CCIndexDoubleVecAttr;

  bool updateCrossSection(const CCStructure &cs, CCIndexDataAttr &indexAttr, const fem::ElasticTriangle3Attr &elementAttr, 
    const CCIndexSet &crossWalls, const IntSet &crackedCells, IntDoubleAttr &crossSectionAttr, 
    IntDoubleAttr &crossSectionWallAttr, CCIndexCCIndexVecAttr &faceList, 
    CCIndexDoubleVecAttr &faceWeightList, CCIndexPoint3dAttr &dirs,CCStructure *csOut = 0, bool verbose = false)
  {
    static const double EPS = 1e-12;

    updateVolumeGeometry(cs, indexAttr);
    if(csOut)
      csOut->clear();

    // First find the cross section directions
    auto &volumes = cs.volumes();
    #pragma omp parallel for mdxNumThreads()
    for(uint i = 0; i < volumes.size(); i++) {
      CCIndex l = volumes[i];
      auto &lDir = dirs[l];
      lDir = Point3d(0,0,0);
      for(CCIndex f : cs.bounds(l)) {
        if(crossWalls.count(f) > 0)
          continue;
        auto &fIdx = indexAttr[f];
        auto &e = elementAttr[f];
        Point3d dir = normalized(fem::triRefCosSinToCurrDir(indexAttr[e.v[0]].pos, indexAttr[e.v[1]].pos, indexAttr[e.v[2]].pos, 
                                                                                                          e.restLen, e.kParCosSin));
        if(dir * lDir < 0)
          dir *= -1;
        lDir += dir * fIdx.measure;
      }
      lDir = normalized(lDir);
    }

//for(CCIndex f : cs.faces())
//  indexAttr[f].selected = false;

    AttrMap<CCIndex, Point3dVec> divideWalls;
    #pragma omp parallel for mdxNumThreads()
    for(uint i = 0; i < volumes.size(); i++) {
      CCIndex vol = volumes[i];
      if(crackedCells.count(indexAttr[vol].label) > 0 and faceList[vol].size() > 0)
        continue;

      //auto &pr = endCenters[vol];
      //Point3d dir = normalized(indexAttr[pr.first].pos - indexAttr[pr.second].pos);
      const auto &dir = dirs[vol];
      SplittingPlane plane(dir, indexAttr[vol].pos);
  
      // neighbouring vertices by vertex
      std::map<CCIndex, std::set<CCIndex>> bdAdj;
    
      // point on edge
      typedef std::pair<CCIndex,double> PointOnEdge;
      std::map<CCIndex, std::vector<PointOnEdge>> faceEdges;
    
      // figure out where to split edges
      for(CCIndex face : cs.incidentCells(vol,2)) {
        CCStructure::CellTuple tuple(cs,face);
        CCIndex firstV = tuple[0];
        std::vector<PointOnEdge> facePt;
        do {
          Point3d p0 = indexAttr[tuple[0]].pos;
          tuple.flip(0);
          Point3d p1 = indexAttr[tuple[0]].pos;
    
          if(plane.inside(p0) != plane.inside(p1)) {
            double t = plane.intersectEdge(p0,p1);
    
            if(t < EPS || t > (1.+EPS)) continue;
            else if(t > (1.-EPS)) {  // line passes through vertex tuple[0]
              facePt.push_back(std::make_pair(tuple[1], (tuple.ro(0) == ccf::POS) ? 0 : 1));
            } else {
              facePt.push_back(std::make_pair(tuple[1], (tuple.ro(0) == ccf::POS) ? (1-t) : t));
            }
          }
          tuple.flip(1);
        } while(tuple[0] != firstV);
    
        // There should be zero, one, or two intersection points in this face.
        switch(facePt.size()) {
          // Zero points means the split doesn't hit the face.
        case 0:
    
          // One means a vertex of the original face,
          // but the split doesn't travel through the face.
        case 1:
          break;
    
          // Two vertices means we have an edge, and may have to split the face.
        case 2:
          // We can check if the new points share an edge;
          // for this to happen, they must both be endpoints!
          {
            CCIndex p0;
            if(facePt[0].second == 0) p0 = cs.edgeBounds(facePt[0].first).first;
            else if(facePt[0].second == 1) p0 = cs.edgeBounds(facePt[0].first).second;
    
            CCIndex p1;
            if(facePt[1].second == 0) p1 = cs.edgeBounds(facePt[1].first).first;
            else if(facePt[1].second == 1) p1 = cs.edgeBounds(facePt[1].first).second;
    
            // If they share an edge, they are adjacent in the membrane boundary.
            if(!(p0.isPseudocell() || p1.isPseudocell()) and cs.dimensionOf(cs.join(p0,p1)) == 1) {
              bdAdj[p0].insert(p1);
              bdAdj[p1].insert(p0);
            }
    
            // Otherwise, they don't make an edge and we mark the face for splitting.
            else {
              if(p0.isPseudocell())
                faceEdges[face].push_back(facePt[0]);
    
              if(p1.isPseudocell())
                faceEdges[face].push_back(facePt[1]);
    
              if(faceEdges[face].empty())
                mdxInfo << QString("splitVolume: two points not an edge??") << endl;
            }
          }
          break;
    
          // Otherwise, something weird has happened.
        default:
          mdxInfo << QString("splitVolume: facePt has %1 vertices").arg(facePt.size()) << ", vol:" << vol << ", face:" << face << endl;
        }
      }
      double area = 0, wallArea = 0;
      auto &divWall = divideWalls[vol];

      auto &fList = faceList[vol];
      auto &fWeightList = faceWeightList[vol]; 
      fList.clear();
      fWeightList.clear();
      for(auto &pr : faceEdges) {
        if(pr.second.size() != 2) {
          mdxInfo << "Bad edge" << endl;
          continue;
        }
//indexAttr[pr.first].selected = true;
        auto &e = elementAttr[pr.first];

        // Find the positions of the cross section segment (deformed)
        CCIndex e1 = (*pr.second.begin()).first;
        double s1 = (*pr.second.begin()).second;
        auto e1b = cs.edgeBounds(e1);
        Point3d p1 = indexAttr[e1b.first].pos * (1.0 - s1) + indexAttr[e1b.second].pos * s1;
  
        CCIndex e2 = (*pr.second.rbegin()).first;
        double s2 = (*pr.second.rbegin()).second;
        auto e2b = cs.edgeBounds(e2);
        Point3d p2 = indexAttr[e2b.first].pos * (1.0 - s2) + indexAttr[e2b.second].pos * s2;

        if(csOut) {
          divWall.push_back(p1);
          divWall.push_back(p2);
        }
  
        // Sum the cell cross section area and wall area
        area += fabs(triangleArea(p1, p2, plane.pointOnPlane));
        double wArea = e.thickness * norm(p1 - p2)/cs.cobounds(pr.first).size();
        wallArea += wArea;

        // Save the faces used
        fList.push_back(pr.first);
        fWeightList.push_back(wArea); // Weight each face by wall area
      }
      int label = indexAttr[vol].label;
      crossSectionAttr[label] = area;
      crossSectionWallAttr[label] = wallArea;
      if(verbose)
        mdxInfo << "Cell:" << label << " has cross section area:" << area << ", wall area:" << wallArea << endl;
    }

    if(csOut) {
      for(auto &pr : divideWalls) {
        CCIndex c = CCIndexFactory.getIndex();
        indexAttr[c].pos = indexAttr[pr.first].pos;
        csOut->addCell(c);
        for(uint i = 0; i < pr.second.size(); i += 2) {
          CCIndex v1 = CCIndexFactory.getIndex();
          CCIndex v2 = CCIndexFactory.getIndex();
          CCIndex e = CCIndexFactory.getIndex();
          csOut->addCell(v1);
          csOut->addCell(v2);
          csOut->addCell(e, +v1-v2);
          indexAttr[v1].pos = pr.second[i];
          indexAttr[v2].pos = pr.second[i+1];

          CCIndex ec = CCIndexFactory.getIndex();
          csOut->addCell(ec, +c-v1);
        }
      }
    }
    return true;
  }

  // Load indices from a file
  CCIndexVec loadCCIndexFromFile(const QString &fileName)
  {
    CCIndexVec loadedCells;

    if(fileName.isEmpty()) {
      mdxInfo << "loadCCIndexFromFile No input file specified" << endl; 
      return loadedCells;
    }

    QFile file(fileName);
    if(!file.open(QIODevice::ReadOnly)) {
      mdxInfo << "loadCCIndexFromFile Cannot open file: " << fileName << endl; 
      return loadedCells;
    }

    QString indices = QString::fromUtf8(file.readAll());
    QStringList indexList = indices.split(QRegExp("\\s+"), QString::SkipEmptyParts);

    for(QString &s : indexList) {
      bool ok = false;
      int idx = s.toInt(&ok);
      if(!ok or idx <= 0)
        continue;
      CCIndex c(idx);
      if(c.isPseudocell())
        continue;
      loadedCells.push_back(c);
    }
    file.close();

    return loadedCells;
  }

  // Separate code from Brendan
  // Helper class DisjointSet (copied from CCFCellStructure.cpp)
  template <typename Index>
  struct DisjointSet
  {
    // DisjointSet type definitions
    typedef std::pair<Index,unsigned int> DisjointData;
    typedef std::map<Index,DisjointData> DisjointMap;
    typedef typename DisjointMap::iterator DisjointIterator;
  
    // DisjointSet member variables
    DisjointMap sets;
  
    // DisjointSet FindSet of iterator
    Index FindSet(DisjointIterator iter)
    {
      if(iter->second.first != iter->first)
        iter->second.first = FindSet(iter->second.first);
      return iter->second.first;
    }
  
    // DisjointSet Link
    void Link(DisjointIterator it1, DisjointIterator it2)
    {
      unsigned int &rank1 = it1->second.second, &rank2 = it2->second.second;
      if(rank1 > rank2) it2->second.first = it1->first;
      else if(rank1 < rank2) it1->second.first = it2->first;
      else { it1->second.first = it2->first; rank2++; }
    }
  
  public:
    // DisjointSet MakeSet
    void MakeSet(Index t)
    { if(sets.count(t) == 0) sets[t] = std::make_pair(t,0); }
  
    // DisjointSet FindSet of Index
    Index FindSet(Index t)
    {
      DisjointIterator iter = sets.find(t);
      if(iter == sets.end()) { MakeSet(t); return t; }
      else return FindSet(iter);
    }
  
    // DisjointSet Union
    void Union(Index t1, Index t2)
    { Link(sets.find(FindSet(t1)), sets.find(FindSet(t2))); }
  
    // DisjointSet Components
    std::set<Index> Components(void)
    {
      std::set<Index> components;
      for(DisjointIterator iter = sets.begin() ; iter != sets.end() ; iter++)
        components.insert(FindSet(iter));
      return components;
    }
  
    // DisjointSet CountComponents
    unsigned int CountComponents(void)
    {
      return Components().size();
    }
  };
  
  bool separateCell(CCStructure &csX, CCIndex cell, const std::vector<CCIndex> &newCells)
  {
    // This only works on CCF cell structures.
    ccf::CCStructure &cs = ensureCCF(csX);
  
    // Utility typedefs.
    typedef CCStructure::FlipI FlipI;
    typedef CCStructure::FlipVectorI FlipVectorI;
  
    // Get cell dimension.
    if(!cs.hasCell(cell)) {
      mdxInfo << QString("separateCell: input cell not in complex") << endl;
      return false;
    }
    CCStructure::Dimension dim = cs.dimensionOf(cell);
  
    // Do nothing with max-dimension cells.
    if(dim == cs.maxDimension()) {
      mdxInfo << "separateCell: can't separate cells of max dimension" << endl;
      return true;
    }
  
    // Get the flips involving the cell.
    FlipVectorI flipsI = cs.matchV(CCIndex::Q, CCIndex::Q, CCIndex::Q, cell);
    FlipVectorI flipsF = cs.matchV(CCIndex::Q, cell, CCIndex::Q, CCIndex::Q);
    FlipVectorI flipsJ = cs.matchV(cell, CCIndex::Q, CCIndex::Q, CCIndex::Q);
  
    // We're going to go through the flips and modify them as needed.
    FlipVectorI newFlips, delFlips;
  
    // (maxdim-1)-cells are handled specially.
    if(dim == cs.maxDimension() - 1) {
      // Consistent INFTY handling requires there to be
      // exactly one cobound flip, with two non-pseudocell cobounds.
      // Otherwise, we just fail.
      if(flipsJ.size() != 1) {
        mdxInfo << "separateCell: separating a cell of max dimension less one "
                << "requires exactly one n-flip, not "
                << flipsJ.size() << endl;
        return false;
      }
  
      if(flipsJ[0].facet[0].isPseudocell() ||
         flipsJ[0].facet[1].isPseudocell()) {
        mdxInfo << "separateCell: cannot separate a cell of max dimension less one "
                << "on the border of the complex" << endl;
        return false;
      }
  
      // Make sure we have exactly one new cell. This one's a fatal error.
      if(newCells.size() != 1) {
        mdxInfo << QString("separateCell: separating cell of dimension %1 " "(max-1) with %2 new cells provided")
                                                                              .arg(dim).arg(newCells.size()) << endl;
        return false;
      }
  
      CCIndex newCell = newCells[0];
  
      // Grab the two facets. Each separated cell is now adjacent to one of them;
      // arbitrarily, cell is adjacent to cob[0], newCell to cob[1].
      CCIndex cob[2] = { flipsJ[0].facet[0], flipsJ[0].facet[1] };
  
      // The boundary of the new cell is the same as the old;
      // accordingly, flips with cell on the interior are just duplicated.
      for(FlipI flip : flipsI) {
        newFlips.push_back(flip.replaceIndex(cell, newCell));
      }
  
      // We only have to redirect the facet-flips through cob[1];
      // in addition, we have to add new flips through INFTY.
      ccf::RO alterRO = -cs.ro(cob[0], CCIndex::TOP);   // precompute
      for(FlipI flip : flipsF) {
        if(flip.interior == cob[1]) {
          delFlips.push_back(flip);
          newFlips.push_back(flip.replaceIndex(cell, newCell));
        }
  
        // Now we handle adding the new INFTY flips.
        bool cob0 = flip.interior == cob[0];
        CCIndex thisCell = cob0 ? cell : newCell;
  
        // We find where cell and newCell flip to by following a rotation
        // until we hit INFTY; if we never do, then they flip to each other.
        CCIndex succ;
        CCStructure::CellTuple tuple(cs, flip.joint, cell, flip.interior);
        do {
          tuple.flip(dim, dim+1);
          if(tuple[dim+1] == CCIndex::INFTY) {
            succ = tuple[dim];
            break;
          }
        } while(tuple[dim] != cell);
        if(succ.isPseudocell())
          succ = cob0 ? newCell : cell;
        else if(cob0) {
          // We're going to be splitting up an existing INFTY flip;
          // we hit this twice (from cob[0] and cob[1]) and delete it only once.
          FlipVectorI vs = cs.matchV(flip.joint, succ, CCIndex::Q, CCIndex::INFTY);
          // There should always be exactly one of these.
          if(vs.size() != 1) {
            mdxInfo << QString("separateCell: should be exactly one INFTY flip " "in near-border INFTY manipulation") << endl;
            return false;
          }
          delFlips.push_back(vs[0]);
        }
  
        // Now, flip is +/- (j, cell, x, cob[i]), and the new flip
        // will be +/- (j, thisCell, succ, INFTY).
        // We replace cell with thisCell, cob[i] with INFTY, and x with succ,
        // but what of the orientations?
        FlipI fInt(flip.replaceIndex(cell, thisCell)
                       .replaceIndex(flip.interior, CCIndex::INFTY)
                       .replaceIndex(flip.otherFacet(cell), succ));
  
        // Note that
        //   1. the cob-flip flipsJ[0] is (cell, cob[0], cob[1], TOP)
        //      so ro(cell, cob[0]) = ro(cob[0], TOP) = -ro(cell, cob[1]),
        //      or ro(cob[0], TOP) = (i==0 ? POS : NEG) * ro(cell, cob[i]);
        //   2. one of the new cob-flips is (cell, cob[0], INFTY, TOP)
        //      and ro(INFTY, TOP) = POS, so ro(cell, INFTY) = NEG,
        //      and similarly ro(newCell, INFTY) = POS;
        //   3. therefore, the target flip is going to be
        //        ro(j, cell) * ro(thisCell, INFTY) * (j, thisCell, succ, INFTY)
        //        = ro(j, cell) * (i==0 ? NEG : POS) * (j, thisCell, succ, INFTY);
        //      the sign of thisCell in this flip is ro(j, cell) * (i==0 ? NEG : POS);
        //   4. but the sign of cell in flip is ro(j, cell) * ro(cell, cob[i])
        //        = ro(j, cell) * (i==0 ? POS : NEG) * ro(cob[0], TOP)
        //   5. setting equal, we use the same orientation as the original exactly when
        //        (i==0 ? NEG : POS) == (i==0 ? POS : NEG) * ro(cob[0], TOP)
        //        ro(cob[0], TOP) == NEG
        //      i.e. multiply by -ro(cob[0], TOP):
        newFlips.push_back(alterRO * fInt);
      }
  
      // The single coboundary flip is removed, and two copies with INFTY
      // in the facets are added.
      {
        FlipI cobFlip = flipsJ[0];
        delFlips.push_back(cobFlip);
        // The first is between cob[0] and INFTY through cell;
        newFlips.push_back(cobFlip.replaceIndex(cob[1], CCIndex::INFTY));
        // the second between INFTY and cob[1] through newCell.
        newFlips.push_back(cobFlip.replaceIndex(cob[0], CCIndex::INFTY)
                                  .replaceIndex(cell, newCell));
      }
    }
  
    // Otherwise, we're dimension-agnostic.
    else {
      // Use our disjoint-set data structure to separate the neighbourhood
      // of the cell and figure out how many parts we'll split into.
      DisjointSet<CCIndex> cobDS;
  
      // An entry for every cell in our coboundary.
      for(const FlipI &flip : flipsF) {
        cobDS.MakeSet(flip.interior);
      }
  
      // Join entries if there's a flip between them.
      for(const FlipI &flip : flipsJ) {
        if(!flip.interior.isPseudocell())
          cobDS.Union(flip.facet[0], flip.facet[1]);
      }
  
      // Get our discrete components.
      std::set<CCIndex> components = cobDS.Components();
      uint numComponents = components.size();
  
      // If there's only one (or zero!) this function doesn't do anything.
      if(numComponents <= 1) {
        mdxInfo << "separateCells: too few components in coboundary of cell" << endl;
        return false;
      }
  
      // We also need as many cell indices as components.
      if(newCells.size() != (numComponents - 1)) {
        mdxInfo << QString("separateCells: cell's coboundary has %1 components, but %2 new cell indices were supplied")
                                          .arg(numComponents).arg(newCells.size()) << endl;
        return false;
      }
  
      // Assign one new cell to each component.
      std::map<CCIndex, CCIndex> cellMap;
      for(CCIndex comp : components) {
        uint n = cellMap.size();
        if(n == 0)
          cellMap[comp] = cell;
        else
          cellMap[comp] = newCells[n-1];
      }
  
      // Replicate the boundary flips.
      for(FlipI flip : flipsI) {
        for(CCIndex newCell : newCells) {
          newFlips.push_back(flip.replaceIndex(cell, newCell));
        }
      }
  
      // Replace entries in the other flips as appropriate.
      for(FlipI flip : flipsF) {
        CCIndex newCell = cellMap[cobDS.FindSet(flip.interior)];
        if(newCell != cell) {
          delFlips.push_back(flip);
          newFlips.push_back(flip.replaceIndex(cell, newCell));
        }
      }
  
      for(FlipI flip : flipsJ) {
        CCIndex newCell = cellMap[cobDS.FindSet(flip.realFacet())];
        if(newCell != cell) {
          delFlips.push_back(flip);
          newFlips.push_back(flip.replaceIndex(cell, newCell));
        }
      }
    }
  
    // If we've made it this far, everything's OK.
    // Add and delete the flips and the new cells and return.
    for(const FlipI &flip : delFlips)
      cs.flips.erase(flip);
  
    cs.flips.insert(newFlips.begin(), newFlips.end());
  
    for(CCIndex ix : newCells)
      cs.dimensionMap.insert(ix, dim);
  
    return true;
  }

  bool SeparateVertices::run(CCStructure &cs, CCIndexDataAttr &indexAttr, const CCIndexVec &sepVs, CCIndexCCIndexUSetAttr &cellMap)
  {
    CCIndexVec sepV, sepE, sepF;
    IntVec countV, countE;

    if(cs.maxDimension() != 3)
      throw QString("%1::run Only works on 3D meshes").arg(name());

    CCIndexSet faces, edges, vertices;

    // Add all the edges and faces attached to the vertices
    for(CCIndex v : sepVs) {
      vertices.insert(v);
      for(CCIndex e : cs.incidentCells(v, 1))
        edges.insert(e);
      for(CCIndex f : cs.incidentCells(v, 2))
        faces.insert(f);
    }
    // Get those to be separated
    for(CCIndex f : faces) {
      if(cs.incidentCells(f, 3).size() < 2)
        continue;
      sepF.push_back(f);
    }
    for(CCIndex e : edges) {
      int lCount = cs.incidentCells(e, 3).size();
      if(lCount < 2)
        continue;
      sepE.push_back(e);
      countE.push_back(lCount - 1);
    }
    for(CCIndex v : vertices) {
      int lCount = cs.incidentCells(v, 3).size();
      if(lCount < 2)
        continue;
      sepV.push_back(v);
      countV.push_back(lCount - 1);
    }

    // Split faces
    int faceC = cs.faces().size();
    for(CCIndex f : sepF) {
      if(cs.cobounds(f).size() != 2) {
        mdxInfo << "Face vols:" << cs.cobounds(f).size() << endl;
        continue;
      }
      CCIndex newf = CCIndexFactory.getIndex();
      if(!separateCell(cs, f, {newf}))
        throw QString("%1::run Failed on face separate").arg(name());
      else {
        cellMap[f].insert(newf);
        indexAttr[newf] = indexAttr[f];
      }
    } 

    // Separate edges
    int edgeC = cs.edges().size();
    for(uint i = 0; i < sepE.size(); i++) {
      CCIndexVec newEs(countE[i]);
      for(uint j = 0; j < newEs.size(); j++)
        newEs[j] = CCIndexFactory.getIndex();
      if(!separateCell(cs, sepE[i], newEs))
        throw QString("%1::run Failed on edge separate").arg(name());
      else
        for(uint j = 0; j < newEs.size(); j++)
          cellMap[sepE[i]].insert(newEs[j]);
    }

    // Separate vertices
    int vertexC = cs.vertices().size();
    for(uint i = 0; i < sepV.size(); i++) {
      CCIndexVec newVs(countV[i]);
      for(uint j = 0; j < newVs.size(); j++)
        newVs[j] = CCIndexFactory.getIndex();
      if(!separateCell(cs, sepV[i], newVs))
        throw QString("%1::run Failed on vertex separate").arg(name());
      else
        for(uint j = 0; j < newVs.size(); j++) {
          cellMap[sepV[i]].insert(newVs[j]);
          indexAttr[newVs[j]] = indexAttr[sepV[i]];
        }
    }

    mdxInfo << cs.faces().size() - faceC << " faces separated" << endl;
    mdxInfo << cs.edges().size() - edgeC << " edges separated" << endl;
    mdxInfo << cs.vertices().size() - vertexC << " vertices separated" << endl;

    return true;
  }

  // Helper function to erase attributes.

  bool FemMembranes::initialize(QWidget *parent)
  {
    mesh = getMesh(parm("Stack"));
    if(!mesh) 
      throw QString("%1::run No stack (mesh)").arg(name());
    ccName = mesh->getCCName(parm("CC Name"));
    if(ccName.isEmpty()) 
      throw QString("%1::run No cell complex").arg(name());
    cs = &mesh->ccStructure(ccName);

    // Get the solver process
    if(parm("Solver Process").isEmpty())
      throw QString("%1::initialize Solver process empty").arg(name());
    if(!getProcess(parm("Solver Process"), solverProcess))
      throw QString("%1::initialize Unable to make solver process: %2").arg(name()).arg(parm("Solver Process"));
    solverProcess->initialize(parent);
    solverStep = solverProcess->parm("Step").toDouble();

    // Get the pressure process if there is one
    if(parm("Pressure Process").isEmpty())
      throw QString("%1::initialize Pressure process empty").arg(name());
    if(!getProcess(parm("Pressure Process"), pressureProcess))
      throw QString("%1::initialize Unable to make pressure process: %2").arg(name()).arg(parm("Pressure Process"));
    pressureProcess->initialize(parent);

    // Get the growth process if there is one
    if(parm("Growth Process").isEmpty() or !getProcess(parm("Growth Process"), growthProcess))
      growthProcess = 0;
    else
      growthProcess->initialize(parent);

    return true;
  }

  bool FemMembranes::step()
  {
    // Grab persistent vars attribute
    if(!mesh)
      throw QString("%1::step No current mesh").arg(name());
    if(!cs)
      throw QString("%1::step No cell complex").arg(name());
    if(!solverProcess)
      throw QString("%1:step Solver process not set").arg(name());
    if(!pressureProcess)
      throw QString("%1:step Pressure process not set").arg(name());

    auto &indexAttr = mesh->indexAttr();
    auto &varsAttr = mesh->attributes().attrMap<QString, QString>("Mesh Vars");
    auto &elementAttr = mesh->attributes().attrMap<CCIndex, fem::ElasticTriangle3>(parm("Element Attribute"));
    auto &materialAttr = mesh->attributes().attrMap<CCIndex, fem::TransIsoMaterial>(parm("Material Attribute"));
    auto &growthAttr = mesh->attributes().attrMap<CCIndex, fem::Growth>(parm("Growth Attribute"));
    double pressure = varsAttr["Pressure"].toDouble();

    // FIXME add serialization for sets to the attributes
    auto &outerWalls = mesh->attributes().attrMap<QString, CCIndexSet>("CCIndex Sets")["Outer Walls"];
    auto &crossWalls = mesh->attributes().attrMap<QString, CCIndexSet>("CCIndex Sets")["Cross Walls"];
    auto &wallType = *mesh->labelMap("Wall Type");
    auto &wallName = mesh->labelName("Wall Type");
    auto &cellType = *mesh->labelMap("Cell Type");

    Point3d longDir = normalized(stringToPoint3d(parm("Long Dir")));
    if(norm(longDir) == 0)
      throw QString("%1::step Bad Long Dir").arg(name());

    auto &crossSectionAttr = mesh->heatAttr<double>("Cross Section", "Labels");
    auto &crossSectionWallAttr = mesh->heatAttr<double>("Cross Section Wall", "Labels");
    auto &faceList = mesh->attributes().attrMap<CCIndex, CCIndexVec>("Face List");
    auto &faceWeightList = mesh->attributes().attrMap<CCIndex, DoubleVec>("Face Weight List");
    auto &cellDirs = mesh->attributes().attrMap<CCIndex, Point3d>("Cell Dirs");

    // On the first step find the crosswalls, set the thickness and calculate the turgor stresses
    if(varsAttr["First Run"].isEmpty()) {
      varsAttr["First Run"] = "No";
      // Mark the end faces
      crossWalls.clear();
      wallName[1] = "Inner Walls";
      wallName[2] = "Outer Walls";
      wallName[3] = "Cross Walls";
      for(auto f : cs->faces()) {
        auto &fIdx = indexAttr[f];
        if(fIdx.label <= 0)
          throw QString("%1::step Faces must be labeled").arg(name());
        // Default type, inner
        wallType[fIdx.label] = 1;

        // Check direction of wall to find cross walls (and ends)
        if(fabs(norm(longDir * indexAttr[f].nrml) - 1.0) < .01) {
          wallType[fIdx.label] = 3;
          crossWalls.insert(f);
        } else if(cs->cobounds(f).size() == 1) {
          // Walls on the outside that are not cross walls are outer walls
          wallType[fIdx.label] = 2;
          outerWalls.insert(f);
        }
      }

      double crossWallThickness = parm("Cross Wall Thickness").toDouble();
      if(crossWallThickness <= 0)
        throw QString("%1::step Cross Wall Thickness must be > 0").arg(name());
      double innerThickness = parm("Inner Thickness").toDouble();
      if(innerThickness < 0)
        throw QString("%1::step Inner thickness must be > 0").arg(name());
      double outerThickness = parm("Outer Thickness").toDouble();
      if(outerThickness < 0)
        throw QString("%1::step Outer thickness must be > 0").arg(name());
      if(outerThickness == 0)
        mdxInfo << "Outer thickenss zero, using existing thickensses" << endl;
      else {
        mdxInfo << "Setting inner thickness to:" << innerThickness << ", outer to:" << outerThickness << endl;
        for(auto &pr : elementAttr)  {
          if(crossWalls.count(pr.first) > 0)
            pr.second.thickness = crossWallThickness;
          else if(outerWalls.count(pr.first) > 0)
            pr.second.thickness = outerThickness;
          else
            pr.second.thickness = innerThickness * 2.0; // Inner walls are shared and doubled
        }
      }

      // Get cross sectional area
      //updateCrossSection(*cs, indexAttr, elementAttr, crossWalls, crossSectionAttr, crossSectionWallAttr, faceList, cellDirs);
    }

    // Run the fem simulation, return true if not converged for another step
    bool converged = !solverProcess->step();

    // Find the cross sections
    QString divideName = parm("Divide CC");
    CCStructure *csDivide = 0;
    if(!divideName.isEmpty()) {
      csDivide = &mesh->ccStructure(divideName);
      mesh->updateAll(divideName);
      csDivide->clear();
    }
    updateCrossSection(*cs, indexAttr, elementAttr, crossWalls, crackedCells, crossSectionAttr, crossSectionWallAttr, faceList, faceWeightList, cellDirs, csDivide, false);

    // Update the stress
    if(parm("Long Stress Signal").isEmpty())
      throw QString("%1::step Dir Stress Signal empty").arg(name());
    auto &longStressAttr = mesh->signalAttr<double>(parm("Long Stress Signal"));
    if(parm("Tissue Stress Signal").isEmpty())
      throw QString("%1::step Tissue Stress Signal empty").arg(name());
    auto &tissueStressAttr = mesh->heatAttr<double>("Tissue Stress", "Labels");

    auto &longStrainAttr = mesh->signalAttr<double>("Long Strain");
    auto &transStrainAttr = mesh->signalAttr<double>("Trans Strain");
    auto &specGrowthAttr = mesh->heatAttr<double>("Spec Growth", "Labels");

    longStressAttr.clear();
    longStrainAttr.clear();
    transStrainAttr.clear();
    tissueStressAttr.clear();
    specGrowthAttr.clear();
    for(const auto &pr : elementAttr) {
      auto &longStress = longStressAttr[pr.first];
      auto &longStrain = longStrainAttr[pr.first];
      auto &transStrain = transStrainAttr[pr.first];
      longStress = longStrain = transStrain = 0;
      if(crossWalls.count(pr.first) > 0) // Skip cross walls
        continue;
 
      auto nrml = normalized(indexAttr[pr.first].nrml);
      auto cb = cs->cobounds(pr.first);
      for(CCIndex l : cb) {
        auto &longDir = cellDirs[l];
        auto transDir = normalized(longDir ^ nrml);
        auto pStress = SymmetricTensor(pr.second.stress).reproject(longDir, transDir, nrml);
        longStress += pStress.evals()[0];
        auto pStrain = SymmetricTensor(pr.second.strain).reproject(longDir, transDir, nrml);
        longStrain += pStrain.evals()[0];
        transStrain += pStrain.evals()[1];
      }
      longStress /= cb.size();
      longStrain /= cb.size();
      transStrain /= cb.size();
    }

    // Find tissue stress
    double outerFactor = parm("Outer Factor").toDouble();
    double innerFactor = parm("Inner Factor").toDouble();
    for(CCIndex l : cs->volumes()) {
      int label = indexAttr[l].label;
      if(label <= 0) {
        mdxInfo << "Bad cell label: " << label << " found for volume:" << l << endl;
        continue;
      }

      // Average wall stress over the cell
      double stress = 0, strain = 0, weight = 0, avgGrowth = 0, avgYoung = 0, avgThresh = 0;
      auto &fList = faceList[l];
      auto &fwList = faceWeightList[l];
      for(uint i = 0; i < fList.size(); i++) {
        CCIndex f = fList[i];
        auto &fIdx = indexAttr[f];

        double wt = fwList[i]; // Weight is the wall area
        stress += wt * longStressAttr[f];
        strain += wt * longStrainAttr[f];

        auto &g = growthAttr[f];
        avgGrowth += wt * g.strainKIso;
        //avgThresh = max(avgThresh, g.threshKIso); // Model uses iso only for growth
        avgThresh += wt * g.threshKIso; // Model uses iso only for growth
        avgYoung += wt * materialAttr[f].youngE1E3; // Soft direction only

        weight += wt;
      }
      if(weight <= 0) {
        mdxInfo << "Bad face weight for volume:" << l << endl;
        continue;
      }
      stress /= weight;
      strain /= weight;
      avgGrowth /= weight;
      avgThresh /= weight;
      avgYoung /= weight;
      double turgorStress = (pressure * crossSectionAttr[label])/crossSectionWallAttr[label];

      if(outerFactor > 0) {
//mdxInfo << "Label:" << label << " turgorStress:" << turgorStress << " stress:" << stress << " inner:" << pressure * innerFactor 
//        << " outer:" << pressure * outerFactor << " factor:" <<  crossSectionAttr[label]/crossSectionWallAttr[label]  << endl;
        if(cellType[label] == 1) {
//mdxInfo << "  inner avgGrowth:" << avgGrowth << " thresh:" << avgThresh << " strain:" << pressure * innerFactor/avgYoung << endl;
          tissueStressAttr[label] = stress - pressure * innerFactor; // Inside cell
          specGrowthAttr[label] = avgGrowth * max(0.0, pressure * innerFactor/avgYoung - avgThresh); // Growth * (strain - thresh)
        } else {
//mdxInfo << "  outer avgGrowth:" << avgGrowth << " thresh:" << avgThresh << " strain:" << pressure * outerFactor/avgYoung << endl;
          if(crackedCells.count(label) > 0)
            tissueStressAttr[label] = 0;
          else
            tissueStressAttr[label] = stress - pressure * outerFactor; // Epidermis
          specGrowthAttr[label] = avgGrowth * max(0.0, pressure * outerFactor/avgYoung - avgThresh); // Growth * (strain - thresh)
        }
      } else {
        if(crackedCells.count(label) > 0)
          tissueStressAttr[label] = 0;
        else
          tissueStressAttr[label] = stress - turgorStress; 
        specGrowthAttr[label] = avgGrowth * max(0.0, turgorStress/avgYoung - avgThresh); // Growth * (strain - thresh)
      }

//mdxInfo << "Spec:" << specGrowthAttr[label] << " avgGrowth:" << avgGrowth << " turgorStress:" << turgorStress << " avgYoung:" 
//                   << avgYoung << " avgThresh:" << avgThresh << " avg strain:" << turgorStress/avgYoung << endl;
    }

    mesh->updateProperties();

    // After cracking we want to step slowly until converged
    if(!converged and !cracking)
      return true;

    // If converged, take a snapshot
    QString snapshotFile = parm("Snapshot File");
    if(!snapshotFile.isEmpty()) {
      int snapshot = varsAttr["Snapshot"].toInt();
      takeSnapshot(QString("%1-%2.jpg").arg(snapshotFile).arg(snapshot++, 5, 10, QChar('0')));
      varsAttr["Snapshot"] = QString("%1").arg(snapshot);
    }
    if(cracking) {
      // After so many crack steps, put the original solver step back
      if(crackSteps-- <= 0 or converged) {
        solverProcess->setParm("Step", QString("%1").arg(solverStep));
        solverProcess->initSolver(cs);
        cracking = false;
        crackDone = true;
      }
      if(!converged)
        return true;
    }

    auto &pressureAttr = mesh->attributes().attrMap<CCIndex, fem::Pressure>(pressureProcess->parm("Pressure Attribute"));

    // First ramp up the pressure
    if(pressure < parm("Pressure Max").toDouble()) {
      pressure += parm("Pressure Step").toDouble();
      pressureProcess->run(*cs, pressureAttr, cs->faces(), pressure, 1);
      mdxInfo << QString("Pressure Step: %1").arg(pressure) << endl;
      varsAttr["Pressure"] = QString("%1").arg(pressure);
      solverProcess->initSolver(cs);

      return true;
    }

    // Finally check for cracks
    double growthTime = varsAttr["Growth Time"].toDouble();
    QString crackVertices = parm("Crack Vertices");
    double crackTime = parm("Crack Time").toDouble();
    if(growthTime >= crackTime and !crackDone and !crackVertices.isEmpty()) {
      cracking = true;
      const auto &vertices = loadCCIndexFromFile(crackVertices);
      if(vertices.size() == 0)
        return false;

      SeparateVertices *separateProcess;
      if(parm("Separate Process").isEmpty() or !getProcess(parm("Separate Process"), separateProcess))
        throw QString("%1::step Unable to make separate process:%2").arg(name()).arg(parm("Separate Process"));

      auto &elementAttr = mesh->attributes().attrMap<CCIndex, fem::ElasticTriangle3>(separateProcess->parm("Element Attribute"));
      auto &materialAttr = mesh->attributes().attrMap<CCIndex, fem::TransIsoMaterial>(separateProcess->parm("Material Attribute"));
      auto &dirichletAttr = mesh->attributes().attrMap<CCIndex, fem::Dirichlet>(separateProcess->parm("Dirichlet Attribute"));
      auto &growthAttr = mesh->attributes().attrMap<CCIndex, fem::Growth>(separateProcess->parm("Growth Attribute"));
      double pressure = separateProcess->parm("Pressure").toDouble();
      int pressureLabel = separateProcess->parm("Pressure Label").toInt();
      separateProcess->run(*cs, mesh->indexAttr(), vertices, 
                    &elementAttr, &materialAttr, &dirichletAttr, &growthAttr, &pressureAttr, pressure, pressureLabel);

      // Change step so we can watch it crack
      solverProcess->setParm("Step", parm("Crack Step"));
      solverProcess->initSolver(cs);
      crackSteps = parm("Crack Steps").toInt();

      mesh->updateAll();

      return true;
    }

    // Get the cracked 3D cell labels if required
    if(crackedCells.size() == 0 and growthTime >= crackTime and !parm("Cracked Cells").isEmpty() and !crackVertices.isEmpty()) {
      for(CCIndex c : loadCCIndexFromFile(parm("Cracked Cells")))
        crackedCells.insert(indexAttr[c].label);
      mdxInfo << "Loading crackedCells:" << crackedCells.size() << endl;
    }

    // Grow 
    if(growthProcess and growthTime < parm("Growth Time").toDouble()) {
      growthProcess->run();
      double growthDt = growthProcess->parm("Growth Dt").toDouble();
      growthTime += growthDt;
      mdxInfo << QString("Growth Step Time %1, step %2").arg(growthTime).arg(growthDt) << endl;
      varsAttr["Growth Time"] = QString("%1").arg(growthTime);

      return true;
    }

    // Otherwise we are done
    return false;
  }

  bool FemMembranes::rewind(QWidget *parent)
  {
    // To rewind, we'll reload the mesh
    Mesh *mesh = currentMesh();
    if(!mesh or mesh->file().isEmpty())
      throw QString("%1::rewind() No current mesh, cannot rewind").arg(name());
    MeshLoad meshLoad(*this);
    meshLoad.setParm("File Name", mesh->file());
    crackDone = false;

    return meshLoad.run();
  }

  REGISTER_PROCESS(FemMembranes);
  REGISTER_PROCESS(FemMembraneSolver);
  REGISTER_PROCESS(FemMembraneDerivs);
  REGISTER_PROCESS(FemMembranePressureDerivs);
  REGISTER_PROCESS(FemMembraneDirichletDerivs);

  REGISTER_PROCESS(FemMembraneRefCfg);
  REGISTER_PROCESS(FemMembraneStressStrain);
  REGISTER_PROCESS(FemMembraneSetMaterial);
  REGISTER_PROCESS(FemMembraneAnisoDir);
  REGISTER_PROCESS(FemMembraneSetPressure);
  REGISTER_PROCESS(FemMembraneSet3DCellPressure);
  REGISTER_PROCESS(FemMembraneSetFacePressureFromVolume);
  REGISTER_PROCESS(FemMembraneSetDirichlet);
  REGISTER_PROCESS(SeparateVertices);

  REGISTER_PROCESS(FemMembraneSetGrowth);
  REGISTER_PROCESS(FemMembraneCellFacesGrowth);
  REGISTER_PROCESS(FemMembraneGrow);

  REGISTER_PROCESS(FemMembraneVisMaterial);
  REGISTER_PROCESS(FemMembraneVisPressure);
  REGISTER_PROCESS(FemMembraneVisGrowth);
  REGISTER_PROCESS(FemMembraneVisDirichlet);
  REGISTER_PROCESS(FemMembraneVisDirections);
  REGISTER_PROCESS(FemMembraneVisThickness);
  REGISTER_PROCESS(FemAnisotropyPropagationFailure);
  REGISTER_PROCESS(SplitEdgesAngle);
  REGISTER_PROCESS(MergeFacesEdgeLength);

  REGISTER_PROCESS(CleanFemAttributes);
  REGISTER_PROCESS(DeleteFemAttributes);

  REGISTER_PROCESS(MakeColorMap);
  REGISTER_PROCESS(AnimateRotation);
}

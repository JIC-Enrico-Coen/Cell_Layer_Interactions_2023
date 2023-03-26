#ifndef FEM_MEMBRANES_HPP
#define FEM_MEMBRANES_HPP

#include <MDXProcessFem.hpp>
#include <MeshProcessSystem.hpp>
#include <MeshProcessStructure.hpp>

namespace mdx
{
  class FemMembraneSolver;
  class FemMembraneGrow;
  class FemMembraneSetPressure;

  class FemMembranes : public Process
  {
  public:
    FemMembranes(const Process &proc) : Process(proc) 
    {
      setName("Model/CCF/00 Fem Membranes");
      setDesc("Fem Simulation of extensometer");

      addParm("Stack", "Stack name, empty for current", "");
      addParm("CC Name", "Cell complex name, empty for current", "");

      addParm("Solver Process", "Name of Fem solver process", "Model/CCF/01 Fem Solver");
      addParm("Growth Process", "Name of the growth process", "Model/CCF/23 Grow");
      addParm("Pressure Process", "Name of the set pressure process", "Model/CCF/14 Set Pressure");
      addParm("Separate Process", "Name of the process to separate vertices", "Model/CCF/18 Separate Vertices");
      addParm("Cross Wall Thickness", "Thickness of cross walls", "5.0");
      addParm("Inner Thickness", "Thickness of inner walls", "0.5");
      addParm("Outer Thickness", "Thickness of outer walls", "10.0");
      addParm("Pressure Max", "Maximum pressure", "0.5");
      addParm("Pressure Step", "Amount to step pressure", "0.025");
      addParm("Growth Time", "Amount of time to grow", "10.0");
      addParm("Crack Time", "Time to crack", "7.0");
      addParm("Crack Vertices", "File name for crack vertices, empty doesn't crack", "");
      addParm("Crack Step", "Time step when cracking", "5.0");
      addParm("Crack Steps", "Number of steps to animate crack", "10");
      addParm("Cracked Cells", "3D cells that get cracked, used to 0 tissue stress", "CrackedCells.txt");
      addParm("Element Attribute", "Name of the element attribute", "Triangle Element");
      addParm("Material Attribute", "Attribute for material", "TransIso Material");
      addParm("Growth Attribute", "Attribute for growth", "Fem Growth");
      addParm("Long Dir", "Longitudinal direction", "0.0 0.0 0.1");
      addParm("Long Stress Signal", "Signal name for longitudinal stress", "Long Stress");
      addParm("Tissue Stress Signal", "Signal name for longitudinal  tissue stress", "Tissue Stress");
      addParm("Divide CC", "Create a cell complex with division walls, empty for none", "Divide Walls");
      addParm("Outer Factor", "Outer cell factor for autonomous stress", "5.4");
      addParm("Inner Factor", "Inner cell factor for autonomous stress", "25.69");

      addParm("Snapshot File", "Snapshot file, empty for no snapshots", "Snapshot/Frame");
    }
    bool initialize(QWidget *parent);
    bool step();
    bool rewind(QWidget *parent);

  private:
    Mesh *mesh = 0;
    QString ccName;
    CCStructure *cs = 0;
    FemMembraneSolver *solverProcess = 0;
    FemMembraneGrow *growthProcess = 0;
    FemMembraneSetPressure *pressureProcess = 0;

    bool crackDone = false;
    bool cracking = false;
    int crackSteps = 0;
    IntSet crackedCells;
    double solverStep = 100.0;
  };

  class FemMembraneSolver : public fem::FemSolver
  {
  public:
    FemMembraneSolver(const Process &proc) : FemSolver(proc) 
    {
      setName("Model/CCF/01 Fem Solver");
      setDesc("Fem solver process");

      // Update parameters with our own defaults
      setParmDefault("Stress-Strain", "Model/CCF/11 Stress-Strain");

      // Add derivatives processes
      addParm("Element Derivs", "Process for element derivatives", "Model/CCF/02 Triangle Derivs");
      addParm("Pressure Derivs", "Process for pressure derivatives", "Model/CCF/03 Pressure Derivs");
      addParm("Dirichlet Derivs", "Process for Dirichlet derivatives", "Model/CCF/04 Dirichlet Derivs");
    }
  };
  class FemMembraneDerivs : public fem::ElementDerivs
  {
  public:
    FemMembraneDerivs(const Process &proc) : ElementDerivs(proc) 
    {
      setName("Model/CCF/02 Triangle Derivs");

      setParmDefault("Element Type", "Linear Triangle");
      setParmDefault("Element Attribute", "Triangle Element");
    }
  };
  class FemMembranePressureDerivs : public fem::PressureDerivs
  {
  public:
    FemMembranePressureDerivs(const Process &proc) : PressureDerivs(proc) 
    {
      setName("Model/CCF/03 Pressure Derivs");
    }
  };
  class FemMembraneDirichletDerivs : public fem::DirichletDerivs
  {
  public:
    FemMembraneDirichletDerivs(const Process &proc) : DirichletDerivs(proc) 
    {
      setName("Model/CCF/04 Dirichlet Derivs");
    }
  };

  class FemMembraneRefCfg : public fem::SetRefCfg
  {
  public:
    FemMembraneRefCfg(const Process &proc) : SetRefCfg(proc) 
    {
      setName("Model/CCF/10 Reference Configuration");

      addParm("Thickness", "Thickness of the membrane elements", "1.0");
      setParmDefault("Element Type", "Linear Triangle");
      setParmDefault("Element Attribute", "Triangle Element");
    }
  };
  class FemMembraneStressStrain : public fem::StressStrain
  {
  public:
    FemMembraneStressStrain(const Process &proc) : StressStrain(proc) 
    {
      setName("Model/CCF/11 Stress-Strain");

      setParmDefault("Element Type", "Linear Triangle");
      setParmDefault("Element Attribute", "Triangle Element");
    }
  };
  class FemMembraneSetMaterial : public fem::SetTransIsoMaterial
  {
  public:
    FemMembraneSetMaterial(const Process &proc) : SetTransIsoMaterial(proc) 
    {
      setName("Model/CCF/12 Set Material Properties");
    }
  };
  class FemMembraneAnisoDir : public fem::SetAnisoDir
  {
  public:
    FemMembraneAnisoDir(const Process &proc) : SetAnisoDir(proc) 
    {
      setName("Model/CCF/13 Set Ansio Dir");

      setParmDefault("Element Type", "Linear Triangle");
      setParmDefault("Element Attribute", "Triangle Element");
    }
  };
  class FemMembraneSetPressure : public fem::SetPressure
  {
  public:
    FemMembraneSetPressure(const Process &proc) : SetPressure(proc) 
    {
      setName("Model/CCF/14 Set Pressure");
    }
  };
  class FemMembraneSet3DCellPressure : public fem::Set3DCellPressure
  {
  public:
    FemMembraneSet3DCellPressure(const Process &proc) : Set3DCellPressure(proc) 
    {
      setName("Model/CCF/15 Set 3D Cell Pressure");
      
      setParmDefault("Pressure Signal", "Cell Pressure");
    }
  };
  class FemMembraneSetFacePressureFromVolume : public fem::SetFacePressureFromVolumes
  {
  public:
    FemMembraneSetFacePressureFromVolume(const Process &proc) : SetFacePressureFromVolumes(proc) 
    {
      setName("Model/CCF/16 Set Face Pressure from Volume");
      
      setParmDefault("Pressure Signal", "Cell Pressure");
    }
  };
  class FemMembraneSetDirichlet : public fem::SetDirichlet
  {
  public:
    FemMembraneSetDirichlet(const Process &proc) : SetDirichlet(proc) 
    {
      setName("Model/CCF/17 Set Dirichlet");
    }
  };

  class FemMembraneSetGrowth : public fem::SetGrowth
  {
  public:
    FemMembraneSetGrowth(const Process &proc) : SetGrowth(proc) 
    {
      setName("Model/CCF/20 Set Growth");
    }
  };
  class FemMembraneCellFacesGrowth : public fem::CellsToFacesGrowth
  {
  public:
    FemMembraneCellFacesGrowth(const Process &proc) : CellsToFacesGrowth(proc) 
    {
      setName("Model/CCF/21 Cells to Faces Growth");
    }
  };
  class FemMembraneGrow : public fem::Grow
  {
  public:
    FemMembraneGrow(const Process &proc) : Grow(proc) 
    {
      setName("Model/CCF/23 Grow");
    }
  };

  class FemMembraneVisMaterial : public fem::VisTransIsoMaterial
  {
  public:
    FemMembraneVisMaterial(const Process &proc) : VisTransIsoMaterial(proc) 
    {
      setName("Model/CCF/30 Visualize Material");
    }
  }; 
  class FemMembraneVisPressure : public fem::VisPressure
  {
  public:
    FemMembraneVisPressure(const Process &proc) : VisPressure(proc) 
    {
      setName("Model/CCF/31 Visualize Pressure");
    }
  };
  class FemMembraneVisGrowth : public fem::VisGrowth
  {
  public:
    FemMembraneVisGrowth(const Process &proc) : VisGrowth(proc) 
    {
      setName("Model/CCF/32 Visualize Growth");
    }
  };
  class FemMembraneVisDirichlet : public fem::VisDirichlet
  {
  public:
    FemMembraneVisDirichlet(const Process &proc) : VisDirichlet(proc) 
    {
      setName("Model/CCF/33 Visualize Dirichlet");
    }
  };
  class FemMembraneVisDirections : public fem::VisDirections
  {
    public:
      FemMembraneVisDirections(const Process &proc) : VisDirections(proc)
      {
        setName("Model/CCF/34 Visualize Directions") ;
      }
  };
  class FemMembraneVisThickness : public Process
  {
  public:
    FemMembraneVisThickness(const Process &proc) : Process(proc)
    {
      setName("Model/CCF/35 Visualize Thickness");
      setDesc("Visualize the element thickness");
      addParm("Stack", "Stack name, empty for current", "");
      addParm("Element Attribute", "Attribute for element", "Triangle Element");
    }

    bool run()
    {
      Mesh *mesh = getMesh(parm("Stack"));
      if(!mesh) 
        throw QString("%1::run No stack (mesh)").arg(name());

      auto &elementAttr = mesh->attributes().attrMap<CCIndex, fem::ElasticTriangle3>(parm("Element Attribute"));
      
      // Add visualization of thickness
      auto &thicknessSignal = mesh->signalAttr<double>("Thickness");
      for(auto &pr : elementAttr)
        thicknessSignal[pr.first] = pr.second.thickness;
      mesh->setSignalBounds(calcBounds(thicknessSignal), "Thickness");
      mesh->setSignal("Thickness");
      mesh->updateProperties();

      return true;
    }
  };

  class FemAnisotropyPropagationFailure : public fem::DisplayFailedAnisotropyPropagation
  { 
    public:
    FemAnisotropyPropagationFailure(const Process &proc) : DisplayFailedAnisotropyPropagation(proc) 
    {
      setName("Model/CCF/36 Display Anisotropy Propagation Failure");
    }
  };

  // Helper class to erase attribute
  template<typename T>
  bool cleanAttr(const CCIndexSet &cells, AttrMap<CCIndex, T> &attrMap)
  {
    CCIndexSet toErase;
    for(auto &pr : attrMap)
      if(cells.count(pr.first) == 0)
        toErase.insert(pr.first);
    for(CCIndex c : toErase)
      attrMap.erase(c);

    return true;
  }

  class CleanFemAttributes : public Process
  {
  public:
    CleanFemAttributes(const Process &proc) : Process(proc) 
    {
      setName("Model/CCF/40 Clean Attributes");
      setDesc("Delete dangling attributes, ie those that do not have a CCIndex in the current mesh");

      addParm("Stack", "Stack name, empty for current", "");
      addParm("CC Name", "Cell complex name, empty for current", "");

      addParm("Element Attribute", "Attribute for element", "Triangle Element");
      addParm("Material Attribute", "Attribute for material", "TransIso Material");
      addParm("Pressure Attribute", "Attribute for pressure", "Fem Pressure");
      addParm("Dirichlet Attribute", "Attribute for Dirichlet conditions", "Fem Dirichlet");
      addParm("Growth Attribute", "Attribute for growth", "Fem Growth");
    }

    bool run()
    {
      Mesh *mesh = getMesh(parm("Stack"));
      if(!mesh) 
        throw QString("%1::run No stack (mesh)").arg(name());
      QString ccName = mesh->getCCName(parm("CC Name"));
      if(ccName.isEmpty()) 
        throw QString("%1::run No cell complex").arg(name());

      // Put all the faces into a set
      auto &cs = mesh->ccStructure(ccName);
      CCIndexSet faces(cs.faces().begin(), cs.faces().end());
      CCIndexSet vertices(cs.vertices().begin(), cs.vertices().end());

      cleanAttr(faces, mesh->attributes().attrMap<CCIndex, fem::ElasticTriangle3>(parm("Element Attribute")));
      cleanAttr(faces, mesh->attributes().attrMap<CCIndex, fem::TransIsoMaterial>(parm("Material Attribute")));
      cleanAttr(faces, mesh->attributes().attrMap<CCIndex, fem::Pressure>(parm("Pressure Attribute")));
      cleanAttr(vertices, mesh->attributes().attrMap<CCIndex, fem::Dirichlet>(parm("Dirichlet Attribute")));
      cleanAttr(faces, mesh->attributes().attrMap<CCIndex, fem::Growth>(parm("Growth Attribute")));

      return true;
    }
  };

  class DeleteFemAttributes : public Process
  {
  public:
    DeleteFemAttributes(const Process &proc) : Process(proc) 
    {
      setName("Model/CCF/41 Delete Attributes");
      setDesc("Delete FEM attributes from the mesh entirely");

      addParm("Stack", "Stack name, empty for current", "");

      addParm("Element Attribute", "Attribute for element", "Triangle Element");
      addParm("Material Attribute", "Attribute for material", "TransIso Material");
      addParm("Pressure Attribute", "Attribute for pressure", "Fem Pressure");
      addParm("Dirichlet Attribute", "Attribute for Dirichlet conditions", "Fem Dirichlet");
      addParm("Growth Attribute", "Attribute for growth", "Fem Growth");
    }

    bool run()
    {
      Mesh *mesh = getMesh(parm("Stack"));
      if(!mesh) 
        throw QString("%1::run No stack (mesh)").arg(name());

      mesh->attributes().erase(parm("Element Attribute"));
      mesh->attributes().erase(parm("Material Attribute"));
      mesh->attributes().erase(parm("Pressure Attribute"));
      mesh->attributes().erase(parm("Dirichlet Attribute"));
      mesh->attributes().erase(parm("Growth Attribute"));

      return true;
    }
  };

  class SeparateVertices : public Process
  {
  public:
    SeparateVertices(const Process& process) : Process(process)
    {
      setName("Model/CCF/18 Separate Vertices");
      setDesc("Separate Vertices");
      setIcon(QIcon(":/images/Sphere.png"));

      addParm("CC Name", "Name of cell complex", "");
      addParm("Element Attribute", "Attribute for element", "Triangle Element");
      addParm("Material Attribute", "Attribute for material", "TransIso Material");
      addParm("Dirichlet Attribute", "Attribute for Dirichlet conditions", "Fem Dirichlet");
      addParm("Growth Attribute", "Attribute for growth", "Fem Growth");
      addParm("Pressure Attribute", "Attribute for pressure", "Fem Pressure");
      addParm("Pressure", "Amount to rest pressure", "0.5");
      addParm("Pressure Label", "Label for pressure", "1");
    }

    bool run()
    {
      Mesh *mesh = currentMesh();
      if(!mesh)
        throw QString("%1::run No mesh selected").arg(name());

      QString ccName = parm("CC Name");
      if(ccName.isEmpty())
        ccName = mesh->ccName();
      if(ccName.isEmpty())
        throw QString("%1::run No CC Name").arg(name());

      CCStructure &cs = mesh->ccStructure(ccName);
      CCIndexDataAttr &indexAttr = mesh->indexAttr();
      auto &elementAttr = mesh->attributes().attrMap<CCIndex, fem::ElasticTriangle3>(parm("Element Attribute"));
      auto &materialAttr = mesh->attributes().attrMap<CCIndex, fem::TransIsoMaterial>(parm("Material Attribute"));
      auto &dirichletAttr = mesh->attributes().attrMap<CCIndex, fem::Dirichlet>(parm("Dirichlet Attribute"));
      auto &growthAttr = mesh->attributes().attrMap<CCIndex, fem::Growth>(parm("Growth Attribute"));
      auto &pressureAttr = mesh->attributes().attrMap<CCIndex, fem::Pressure>(parm("Pressure Attribute"));
      double pressure = parm("Pressure").toDouble();
      int pressureLabel = parm("Pressure Label").toInt();
      mesh->updateAll();
      return run(cs, indexAttr, selectedVertices(cs, indexAttr), 
                       &elementAttr, &materialAttr, &dirichletAttr, &growthAttr, &pressureAttr, pressure, pressureLabel);
    }

    bool run(CCStructure &cs, CCIndexDataAttr &indexAttr, const CCIndexVec &vertices, fem::ElasticTriangle3Attr *elementAttr = 0, 
              fem::TransIsoMaterialAttr *materialAttr = 0, fem::DirichletAttr *dirichletAttr = 0, fem::GrowthAttr *growthAttr = 0, 
              fem::PressureAttr *pressureAttr = 0, double pressure = 0, int pressureLabel = 1)
    {
      CCIndexCCIndexUSetAttr cellMap;
      if(!run(cs, indexAttr, vertices, cellMap))
        return false;

      // Update Attributes

      // Elements
      if(elementAttr) {
        // FIXME Can we do this some other way? Need something like a subdivision object
        Mesh *mesh = currentMesh();
        auto &crossWalls = mesh->attributes().attrMap<QString, CCIndexSet>("CCIndex Sets")["Cross Walls"];
        
        for(auto &pr : *elementAttr) {
          auto &cells = cellMap[pr.first];
          if(cells.size() == 0)
            continue;
          if(cells.size() > 1) {
            mdxInfo << "Face split into more than 2 faces" << endl;
            continue;
          }
          CCIndex f = *cells.begin();
          auto &fEle = (*elementAttr)[f];

          // Copy the element data
          fEle = pr.second;

          // Split the thickness
          fEle.thickness /= 2.0;
          pr.second.thickness /= 2.0;

          // Update cross wall if it is one
          if(crossWalls.count(pr.first) > 0)
            crossWalls.insert(f);
        }

        // Map new vertices to old
        CCIndexCCIndexAttr vMap;
        for(CCIndex v : cs.vertices()) {
          auto it = cellMap.find(v);
          if(it != cellMap.end())
            for(CCIndex n : it->second)
              vMap[n] = v;
        }

        // Update any incorrect vertices
        for(auto &pr : *elementAttr) {
          CCIndexCCIndexAttr fMap;
          for(CCIndex v : cs.faceVertices(pr.first)) {
            auto it = vMap.find(v);
            if(it == vMap.end())
              fMap[v] = v;
            else
              fMap[vMap[v]] = v;
          }
          for(int i = 0; i < 3; i++)
            pr.second.v[i] = fMap[pr.second.v[i]];
        }
      }

      if(materialAttr) {
        for(auto &pr : *materialAttr)
          for(auto &c : cellMap[pr.first])
            (*materialAttr)[c] = pr.second;
      }

      if(dirichletAttr) {
        for(auto &pr : *dirichletAttr)
          for(auto &c : cellMap[pr.first])
            (*dirichletAttr)[c] = pr.second;
      }

      if(growthAttr) {
        for(auto &pr : *growthAttr)
          for(auto &c : cellMap[pr.first])
            (*growthAttr)[c] = pr.second;
      }

      if(pressureAttr) {
        // Pressure needs to be reset for all the faces
        FemMembraneSetPressure(*this).run(cs, *pressureAttr, cs.faces(), pressure, pressureLabel); // FIXME Should do only changed faces?
      }

      return true;
    }
    bool run(CCStructure &cs, CCIndexDataAttr &indexAttr, const CCIndexVec &vertices, CCIndexCCIndexUSetAttr &cellMap);
  };

  class SplitEdgesAngle : public SplitEdges
  {
  public:
    SplitEdgesAngle(const Process &proc) : SplitEdges(proc) 
    {
      setName("Model/CCF/37 Split Edges by Angle");
      setDesc("Split active edges within an angle of a direction");

      addParm("Stack", "Stack name, empty for current", "");
      addParm("CC Name", "Cell complex name, empty for current", "");
      addParm("Direction", "Direction of edges", "0 0 1");
      addParm("Angle", "Cosine of tolerence angle", ".9");
    }

    bool run()
    {
      Mesh *mesh = getMesh(parm("Stack"));
      if(!mesh) 
        throw QString("%1::run No stack (mesh)").arg(name());
      QString ccName = mesh->getCCName(parm("CC Name"));
      if(ccName.isEmpty()) 
        throw QString("%1::run No cell complex").arg(name());

      Point3d dir = normalized(stringToPoint3d(parm("Direction")));
      double cosAngle = parm("Angle").toDouble();

      auto &cs = mesh->ccStructure(ccName);
      auto &indexAttr = mesh->indexAttr();

      CCIndexVec edges;
      for(CCIndex e : activeEdges(cs, indexAttr)) {
        auto eb = cs.edgeBounds(e);
        Point3d eDir = indexAttr[eb.first].pos - indexAttr[eb.second].pos; 
        if(fabs(dir * normalized(eDir)) >= cosAngle)
          edges.push_back(e);
      }

      mesh->updateAll();
      return SplitEdges::run(cs, indexAttr, edges);
    }
  };

  class MakeColorMap : public Process
  {
  public:
    MakeColorMap(const Process &proc) : Process(proc) 
    {
      setName("Model/CCF/99 Make Color Map");
      setDesc("Convert list of float colors to color map");

      addParm("Input File", "Name of input file", "input.txt");
      addParm("Output File", "Name of input file", "output.mdxc");
    }

    bool run()
    {
      QFile infile(parm("Input File"));
      if(!infile.open(QIODevice::ReadOnly))
        throw QString("%1 Cannot open file").arg(name());

      DoubleVec colors;
      QString indices = QString::fromUtf8(infile.readAll());
      for(auto s : indices.split(QRegExp("\\s+"), QString::SkipEmptyParts))
        colors.push_back(s.toDouble());
      int sz = colors.size()/3;
      infile.close();

      QFile outfile(parm("Output File"));
      if(!outfile.open(QIODevice::WriteOnly))
        throw QString("%1 Cannot open file").arg(name());
      QTextStream out(&outfile);
      out << "MDX COLORS 1.0" << endl;
      out << sz << endl;
      for(uint i = 0; i < sz; i++) {
        int j = i * 3;
        out << QString("#%1%2%3").arg(std::min(uint(colors[j]*256),255u), 2, 16, QChar('0'))
                                 .arg(std::min(uint(colors[j+1]*256),255u), 2, 16, QChar('0'))
                                 .arg(std::min(uint(colors[j+2]*256), 255u), 2, 16, QChar('0')) << endl;
      }
      outfile.close();

      return true;
    }
  };

  class AnimateRotation : public Process
  {
  public:
    AnimateRotation(const Process &proc) : Process(proc) 
    {
      setName("Model/CCF/98 Make Animation");
      setDesc("Animate movement");

      addParm("Stack", "Stack name, empty for current", "");
      addParm("CC Name", "Cell complex name, empty for current", "");
      addParm("Rotation", "Rotation in degrees", "360");
      addParm("Vector", "Vector to rotate around", "0 0 1");
      addParm("Steps", "Steps to perform the rotation", "100");
      addParm("Snapshot File", "Snapshot file, empty for no snapshots", "Snapshot/Frame");
      addParm("Start Frame", "Starting number for frame", "0");
      addParm("Translate Origin", "Translate to center of mass", "Yes", booleanChoice());
    }

    bool initialize(QWidget *parent)
    {
      Stack *stack = getStack(parm("Stack"));
      if(!stack) 
        throw QString("%1::run No stack").arg(name());

      Mesh *mesh = getMesh(stack->name());
      if(!mesh) 
        throw QString("%1::run No mmesh").arg(name());
      QString ccName = mesh->getCCName(parm("CC Name"));
      if(ccName.isEmpty()) 
        throw QString("%1::run No cell complex").arg(name());

      auto &cs = mesh->ccStructure(ccName);
      auto &indexAttr = mesh->indexAttr();

      double rotation = parm("Rotation").toDouble();
      int steps = parm("Steps").toDouble();
      if(steps < 1)
        throw QString("%1::initialize Steps must be >= 1").arg(name());

      int startFrame = parm("Start Frame").toInt();
      int pauseFrames = parm("Pause Frames").toInt();

      QString snapshotFile = parm("Snapshot File");
      auto &frame = stack->frame();
      double a = (rotation/180.0) * M_PI/double(steps);
      Point3d v = normalized(stringToPoint3d(parm("Vector")));
      if(fabs(norm(v) - 1) > .00001)
        throw QString("%1::initialize Bad vector %2").arg(name()).arg(parm("Vector"));

      // Define the rotation, and translation if required 
      const Quaternion qInit(v, a);
      qglviewer::Quaternion q(qInit.x(), qInit.y(), qInit.z(), qInit.w());
      qglviewer::Vec c;

      bool translate = stringToBool(parm("Translate Origin"));
      double weight = 0;
      CCIndexDataAttr tmpIndexAttr;
      auto tmpPosition = frame.position();
      if(translate) {
        Point3d center;
        for(CCIndex f : cs.faces()) {
          auto &fIdx = indexAttr[f];
          center += fIdx.pos * fIdx.measure;
          weight += fIdx.measure;
        }
        if(weight <= 0)
          throw QString("%1::initialize No faces in mesh for translation").arg(name());
        center /= weight;

        // Set the translation
        c = frame.inverseCoordinatesOf(qglviewer::Vec(center));
      }
      auto position = frame.position();
      for(int i = 0; i <= steps; i++) {
        QCoreApplication::processEvents();
        if(!progressAdvance())
          return false;
        if(translate)
          frame.rotateAroundPoint(q, c);
        else
          frame.rotate(q);
        updateState();
        updateViewer();
        if(!snapshotFile.isEmpty())
          takeSnapshot(QString("%1-%2.jpg").arg(snapshotFile).arg(startFrame++, 5, 10, QChar('0')));
      }

      return true;
    }
  };

  class MergeFacesEdgeLength : public Process
  {
  public:
    MergeFacesEdgeLength(const Process &proc) : Process(proc) 
    {
      setName("Model/CCF/38 Merge Faces Edge Length");
      setDesc("Merge faces by edge length");

      addParm("Stack", "Stack name, empty for current", "");
      addParm("CC Name", "Cell complex name, empty for current", "");
      addParm("Length", "Length of edges to merge", "8.1");
    }

    bool run()
    {
      Mesh *mesh = getMesh(parm("Stack"));
      if(!mesh) 
        throw QString("%1::run No stack (mesh)").arg(name());
      QString ccName = mesh->getCCName(parm("CC Name"));
      if(ccName.isEmpty()) 
        throw QString("%1::run No cell complex").arg(name());

      double length = parm("Length").toDouble();

      auto &cs = mesh->ccStructure(ccName);
      auto &indexAttr = mesh->indexAttr();

      mesh->updateAll();

      CCIndexVec edges;
      for(CCIndex e : activeEdges(cs, indexAttr)) {
        auto eb = cs.edgeBounds(e);
        if(norm(indexAttr[eb.first].pos - indexAttr[eb.second].pos) >= length) {
          edges.push_back(e);
          indexAttr[eb.first].selected = indexAttr[eb.second].selected = true;
        }
      }

      for(CCIndex e : edges) {
        CCIndexVec faces;
        for(CCIndex f : cs.cobounds(e))
          faces.push_back(f);
        if(faces.size() == 2)
          mergeFaces(cs, indexAttr, faces);
      }

      return true;
    }
  };

//  class Separate3DCellsSelect : public Process
//  {
//  public:
//    Separate3DCellsSelect(const Process& process) : Process(process)
//    {
//      setName("Model/CCF/18 Separate 3D Cells");
//      setDesc("Separate 3D cells by duplicating vertices and faces");
//      setIcon(QIcon(":/images/Sphere.png"));
//
//      addParm("CC Name", "Name of cell complex", "");
//      addParm("CC Out", "Name of output cell complex", "");
//      addParm("2D Mesh", "Make a 2D mesh", "No", booleanChoice());
//      addParm("Use Face Selection", "Only separate selected faces", "Yes", booleanChoice());
//
//      addParm("Element Attribute", "Attribute for element", "Triangle Element");
//      addParm("Material Attribute", "Attribute for material", "TransIso Material");
//      addParm("Dirichlet Attribute", "Attribute for Dirichlet conditions", "Fem Dirichlet");
//      addParm("Pressure Attribute", "Attribute for pressure", "Fem Pressure");;
//      addParm("Growth Attribute", "Attribute for growth", "Fem Growth");;
//    }
//
//    bool run()
//    {
//      Mesh *mesh = currentMesh();
//      if(!mesh)
//        throw QString("%1::run No mesh selected").arg(name());
//
//      QString ccName = parm("CC Name");
//      if(ccName.isEmpty())
//        ccName = mesh->ccName();
//      if(ccName.isEmpty())
//        throw QString("%1::run No CC Name").arg(name());
//
//      QString ccOut = parm("CC Out");
//      if(ccOut.isEmpty())
//        ccOut = ccName + " Separated";
//
//      CCStructure &cs = mesh->ccStructure(ccName);
//      CCStructure &csOut = mesh->ccStructure(ccOut);
//      CCIndexDataAttr &indexAttr = mesh->indexAttr();
//      bool mesh2D = stringToBool(parm("2D Mesh"));
//      csOut = mesh2D ? CCStructure(2) : CCStructure(3);
//      bool useSelect = stringToBool(parm("Use Face Selection"));
//      if(useSelect and mesh2D)
//        throw QString("%1::run All faces need to be separated to make a 2D mesh (not just selected)").arg(name());
//
//      CCIndexCCIndexUSetAttr cellMap;
//      bool result = run(cs, csOut, indexAttr, cellMap, useSelect);
//      if(result) {
//        // Update Attributes
//        auto &elementAttr = mesh->attributes().attrMap<CCIndex, fem::ElasticTriangle3>(parm("Element Attribute"));
//        CCIndexUSet updated;
//        for(auto &pr : elementAttr) {
//          for(auto &f : cellMap[pr.first]) {
//            auto fVols = csOut.incidentCells(f, 3);
//            elementAttr[f].setVertices(csOut, indexAttr, f);
//            //elementAttr[f].setRefCfgToCurrent(indexAttr, true);
//            elementAttr[f].thickness = pr.second.thickness;
//            elementAttr[f].restLen = pr.second.restLen;
//          }
//        }
//
//        auto &materialAttr = mesh->attributes().attrMap<CCIndex, fem::TransIsoMaterial>(parm("Material Attribute"));
//        for(auto &pr : materialAttr)
//          for(auto &c : cellMap[pr.first])
//            materialAttr[c] = pr.second;
//
//        auto &dirichletAttr = mesh->attributes().attrMap<CCIndex, fem::Dirichlet>(parm("Dirichlet Attribute"));
//        for(auto &pr : dirichletAttr)
//          for(auto &c : cellMap[pr.first])
//            dirichletAttr[c] = pr.second;
//
//        auto &growthAttr = mesh->attributes().attrMap<CCIndex, fem::Growth>(parm("Growth Attribute"));
//        for(auto &pr : growthAttr)
//          for(auto &c : cellMap[pr.first])
//            growthAttr[c] = pr.second;
//
//        // Now do pressure for the whole mesh
//        auto &pressureAttr = mesh->attributes().attrMap<CCIndex, fem::Pressure>(parm("Pressure Attribute"));
//        FemMembraneSetPressure(*this).run(csOut, pressureAttr, csOut.faces(), 0.5, 1); // FIXME Does not use pressure in Gui 
//
//        cs = csOut;
//        mesh->erase(ccOut);
//        mesh->setCCName(ccName);
//        CleanFemAttributes(*this).runWithCurrentParms();
//            
//      }
//
//      return result;
//    }
//    bool run(const CCStructure &cs, CCStructure &csOut, CCIndexDataAttr &indexAttr, CCIndexCCIndexUSetAttr &cellMap, bool useSelect = true);
//  };
}
#endif


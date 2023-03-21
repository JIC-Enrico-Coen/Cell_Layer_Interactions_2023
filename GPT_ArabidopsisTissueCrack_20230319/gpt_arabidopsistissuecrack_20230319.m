function [m,result] = gpt_arabidopsistissuecrack_20230319( m, varargin )
%[m,result] = gpt_arabidopsistissuecrack_20230319( m, varargin )
%   Morphogen interaction function.
%   Written at 2023-03-20 10:33:12.
%   GFtbox revision 20211011, 2020-10-11 11:00.

% The user may edit any part of this function lying between lines that
% begin "%%% USER CODE" and "%%% END OF USER CODE".  Those lines themselves
% delimiters themselves must not be moved, edited, deleted, or added.

    result = [];
    if isempty(m), return; end

    setGlobals();
    
    % Handle new-style callbacks.
    if nargin > 1
        if exist('ifCallbackHandler','file')==2
            [m,result] = ifCallbackHandler( m, varargin{:} );
        end
        return;
    end

    fprintf( 1, '%s found in %s\n', mfilename(), which(mfilename()) );

    realtime = m.globalDynamicProps.currenttime;
    dt = m.globalProps.timestep;

%%% USER CODE: INITIALISATION
    if (Steps(m)==0) && m.globalDynamicProps.doinit
        % Put any code here that should only be performed at the start of
        % the simulation.

        % Reset several fields of m to their default states.
        % Give the command "help resetMeshValues" for details.
        % You can delete this if you do not want this to happen.
        m = resetMeshValues( m );
        
        m.userdata.numelementrings = 6;
        m.userdata.numvxtiers = 13;
        m.userdata.zmin = min(m.FEnodes(:,3));
        m.userdata.zmax = max(m.FEnodes(:,3));
        m.userdata.height = m.userdata.zmax - m.userdata.zmin;
        m.userdata.vxringheights = (0:(m.userdata.numvxtiers-1))*(m.userdata.height/(m.userdata.numvxtiers-1)) - (m.userdata.numvxtiers-1)/2;
        m.userdata.radii = sqrt( sum( m.FEnodes(:,[1 2]).^2, 2 ) );
        m.userdata.radius = max( m.userdata.radii );
        m.userdata.angles = atan2( m.FEnodes(:,2), m.FEnodes(:,1) );
        m.userdata.numsplits = [];
        m.outputcolors.residualgrowth = [ 0 0 1; 1 0 0 ]; % Red for tension, blue for compression.
    end
    m = setProjectOptions( m );
    OPTIONS = getModelOptions( m );
    printModelOptions( m );
    
    xxxx = 1;
%%% END OF USER CODE: INITIALISATION

%%% SECTION 1: ACCESSING MORPHOGENS AND TIME.
%%% AUTOMATICALLY GENERATED CODE: DO NOT EDIT.

% Each call of getMgenLevels below returns four results:
% XXX_i is the index of the morphogen called XXX.
% XXX_p is the vector of all of its values.
% XXX_a is its mutation level.
% XXX_l is the "effective" level of the morphogen, i.e. XXX_p*XXX_a.
% In SECTION 3 of the automatically generated code, all of the XXX_p values
% will be copied back into the mesh.

    [kpar_i,kpar_p,kpar_a,kpar_l] = getMgenLevels( m, 'KPAR' );  %#ok<ASGLU>
    [kpar2_i,kpar2_p,kpar2_a,kpar2_l] = getMgenLevels( m, 'KPAR2' );  %#ok<ASGLU>
    [kper_i,kper_p,kper_a,kper_l] = getMgenLevels( m, 'KPER' );  %#ok<ASGLU>
    [pol_i,pol_p,pol_a,pol_l] = getMgenLevels( m, 'POL' );  %#ok<ASGLU>
    [pol2_i,pol2_p,pol2_a,pol2_l] = getMgenLevels( m, 'POL2' );  %#ok<ASGLU>
    [id_epidermis_i,id_epidermis_p,id_epidermis_a,id_epidermis_l] = getMgenLevels( m, 'ID_EPIDERMIS' );  %#ok<ASGLU>
    [id_inner_i,id_inner_p,id_inner_a,id_inner_l] = getMgenLevels( m, 'ID_INNER' );  %#ok<ASGLU>
    [v_tension_i,v_tension_p,v_tension_a,v_tension_l] = getMgenLevels( m, 'V_TENSION' );  %#ok<ASGLU>
    [s_strainret_i,s_strainret_p,s_strainret_a,s_strainret_l] = getMgenLevels( m, 'S_STRAINRET' );  %#ok<ASGLU>
    [id_split_i,id_split_p,id_split_a,id_split_l] = getMgenLevels( m, 'ID_SPLIT' );  %#ok<ASGLU>
    [v_ring_i,v_ring_p,v_ring_a,v_ring_l] = getMgenLevels( m, 'V_RING' );  %#ok<ASGLU>
    [v_level_i,v_level_p,v_level_a,v_level_l] = getMgenLevels( m, 'V_LEVEL' );  %#ok<ASGLU>
    [v_fracture_i,v_fracture_p,v_fracture_a,v_fracture_l] = getMgenLevels( m, 'V_FRACTURE' );  %#ok<ASGLU>
    [s_weakness_i,s_weakness_p,s_weakness_a,s_weakness_l] = getMgenLevels( m, 'S_WEAKNESS' );  %#ok<ASGLU>
    [id_surface_i,id_surface_p,id_surface_a,id_surface_l] = getMgenLevels( m, 'ID_SURFACE' );  %#ok<ASGLU>
    [id_interior_i,id_interior_p,id_interior_a,id_interior_l] = getMgenLevels( m, 'ID_INTERIOR' );  %#ok<ASGLU>

% Mesh type: volumetric
%          FEtype: 84.00000081
%        axisdivs: 12
%      circumdivs: 36
%       coneangle: 0
%         dealign: 0
%          hollow: 0
%       innerdivs: 0
%             new: 1
%        position: 0
%           rings: 6
%            size: 5.00000012

%            Morphogen    Diffusion   Decay   Dilution   Mutant
%            --------------------------------------------------
%                 KPAR         ----    ----       ----     ----
%                KPAR2         ----    ----       ----     ----
%                 KPER         ----    ----       ----     ----
%                  POL         ----    ----       ----     ----
%                 POL2         ----    ----       ----     ----
%         ID_EPIDERMIS         ----    ----       ----     ----
%             ID_INNER         ----    ----       ----     ----
%            V_TENSION         ----    ----       ----     ----
%          S_STRAINRET         ----    ----       ----     ----
%             ID_SPLIT         ----    ----       ----     ----
%               V_RING         ----    ----       ----     ----
%              V_LEVEL         ----    ----       ----     ----
%           V_FRACTURE         ----    ----       ----     ----
%           S_WEAKNESS         ----    ----       ----     ----
%           ID_SURFACE         ----    ----       ----     ----
%          ID_INTERIOR         ----    ----       ----     ----


%%% USER CODE: MORPHOGEN INTERACTIONS

% In this section you may modify the mesh in any way that does not
% alter the set of nodes.

    if (Steps(m)==0) && m.globalDynamicProps.doinit
        numrings = 6;
        outersteps = numrings*6;
        radius = max( max( m.FEnodes(:,[1 2]) ) );
        
        id_surface_p = double( m.userdata.radii >= radius * (cos(pi/outersteps) - 0.0001) );
        id_interior_p = 1 - id_surface_p;
%         id_epidermis_pX = double( m.userdata.radii > radius * OPTIONS.coreradius );
        id_epidermis_p = double( m.userdata.radii > radius * (OPTIONS.corerings/numrings) * (cos(pi/(OPTIONS.corerings*6)) * (1-0.0001)) );
        id_inner_p = 1 - id_epidermis_p;
        pol_p = max(m.FEnodes(:,3)) - m.FEnodes(:,3);
        
        if OPTIONS.outerbulkmodulus ~= 1
            outerElements = perVertextoperFE( m, id_epidermis_p, 'min' )==1;
            m.cellbulkmodulus( outerElements ) = OPTIONS.outerbulkmodulus;
            m.cellstiffness = IsotropicStiffnessMatrix( m.cellbulkmodulus, m.cellpoisson, m.globalProps.plasticGrowth );
            
            % Check.
%             outerVxs = double( perFEtoperVertex( m, double(outerElements) ) > 0 );
%             all( outerVxs==id_epidermis_p )
        end
        
        m = leaf_setmorphogenrole( m, 'STRAINRET', s_strainret_i );
        
        % V_LEVEL labels all vertexes on the horizontal planes of vertexes
        % with their level, from 1 (the bottom surface) to 13 (the top
        % surface). Vertexes not on these planes are given value 0.
        v_level_p = (m.FEnodes(:,3) - m.userdata.zmin) * ((m.userdata.numvxtiers-1)/m.userdata.height) + 1;
        v_level_p( abs(v_level_p - round(v_level_p)) > 0.01 ) = 0;
        v_level_p = round(v_level_p);
        
        v_ring_p = m.userdata.radii * (m.userdata.numelementrings/m.userdata.radius);
        v_ring_p( abs(v_ring_p - round(v_ring_p)) > 0.001 ) = 0;
        v_ring_p = round(v_ring_p);
        
        s_weakness_p = rand( getNumberOfVertexes(m), 1 );
        m = leaf_mgen_conductivity( m, s_weakness_i, 0.1 );
    end
    
    m = leaf_mgen_plotpriority( m, {'id_split','id_epidermis'}, [1 1]  );

    if meshAtOrBeforeTime( m, OPTIONS.diffusiontime )
        % Rescale s_weakness_p to the range 0 ... OPTIONS.weakness.
        s_weakness_p = s_weakness_p - min(s_weakness_p);
        s_weakness_p = s_weakness_p * OPTIONS.weakness/max(s_weakness_p);
    end
    
    if meshAtTime( m, OPTIONS.diffusiontime )
        % Turn off diffusion for S_WEAKNESS.
        m = leaf_mgen_conductivity( m, s_weakness_i, 0 );
        
        % Start growth, for KPAR only. KPAR2 and KPER are always zero.
        kpar_p = OPTIONS.innergrowthrate * id_inner_p + OPTIONS.outergrowthrate * (1 - id_inner_p);
        kallperp = OPTIONS.innerhorizgrowthrate * id_inner_p + OPTIONS.outerhorizgrowthrate * (1 - id_inner_p);
        kpar2_p = kallperp;
        kper_p = kallperp;
    end
    
    s_strainret_p(:) = OPTIONS.strainretention;
%%% END OF USER CODE: MORPHOGEN INTERACTIONS

%%% SECTION 3: INSTALLING MODIFIED VALUES BACK INTO MESH STRUCTURE
%%% AUTOMATICALLY GENERATED CODE: DO NOT EDIT.

    m.morphogens(:,kpar_i) = kpar_p;
    m.morphogens(:,kpar2_i) = kpar2_p;
    m.morphogens(:,kper_i) = kper_p;
    m.morphogens(:,pol_i) = pol_p;
    m.morphogens(:,pol2_i) = pol2_p;
    m.morphogens(:,id_epidermis_i) = id_epidermis_p;
    m.morphogens(:,id_inner_i) = id_inner_p;
    m.morphogens(:,v_tension_i) = v_tension_p;
    m.morphogens(:,s_strainret_i) = s_strainret_p;
    m.morphogens(:,id_split_i) = id_split_p;
    m.morphogens(:,v_ring_i) = v_ring_p;
    m.morphogens(:,v_level_i) = v_level_p;
    m.morphogens(:,v_fracture_i) = v_fracture_p;
    m.morphogens(:,s_weakness_i) = s_weakness_p;
    m.morphogens(:,id_surface_i) = id_surface_p;
    m.morphogens(:,id_interior_i) = id_interior_p;

%%% USER CODE: FINALISATION

% In this section you may modify the mesh in any way whatsoever.
        
        
    if (Steps(m)==0) && m.globalDynamicProps.doinit && OPTIONS.makeinitcracks
        % Create some initial cracks.
        
        if true
            % Make random cracks.
            margin = 2;
            m.userdata.crackvxs = minimalRandomCracks( m, m.userdata.numvxtiers, m.userdata.numvxtiers - 2*margin - 3, margin, 2 );
%         elseif true
%             % Make random cracks.
%             margin = 2;
%             m.userdata.crackvxs = findRandomCracks( m, m.userdata.numvxtiers, m.userdata.numvxtiers - 2*margin - 3, margin, 2 );
        else
            % Make a single crack.
            radiusOk = (m.userdata.radii > OPTIONS.crackradius);
            heightOk = abs(m.FEnodes(:,3) - OPTIONS.crackheight) < 0.0001;
            angle = atan2( m.FEnodes(:,2), m.FEnodes(:,1) );
            angleOk = (angle > OPTIONS.crackstartangle) & (angle < OPTIONS.crackendangle);
            m.userdata.crackvxs = radiusOk & heightOk & angleOk;
        end
        
        id_split_i = FindMorphogenIndex2( m, 'id_split' );
        m.morphogens(:,id_split_i) = double( m.userdata.crackvxs );
        m = makeCrack( m, m.userdata.crackvxs );
    end
%%% END OF USER CODE: FINALISATION

end

function [m,result] = ifCallbackHandler( m, fn, varargin )
    result = [];
    if exist(fn,'file') ~= 2
        return;
    end
    [m,result] = feval( fn, m, varargin{:} );
end


%%% USER CODE: SUBFUNCTIONS

% Here you may write any functions of your own, that you want to call from
% the interaction function, but never need to call from outside it.
% Remember that they do not have access to any variables except those
% that you pass as parameters, and cannot change anything except by
% returning new values as results.
% Whichever section they are called from, they must respect the same
% restrictions on what modifications they are allowed to make to the mesh.

% The GFtbox_..._Callback routines can be deleted if you do not use them.
% Those that you retain will be automatically called by GFtbox at certain
% points in the simulation cycle.
% If you retain them, their headers specifying their arguments and results
% must not be altered.

function [m,result] = GFtbox_Precelldivision_Callback( m, ci ) %#ok<DEFNU>
    result = [];
    % Your code here.

    % If a nonempty result is to be returned, it should be a struct
    % with fields result.divide, result.dividepoint, and result.perpendicular.
end

function [m,result] = GFtbox_Postcelldivision_Callback( m, ci, cei, newci, newcei, oe1, oe2, ne1, ne2, ne3 ) %#ok<DEFNU>
    result = [];
    % Your code here.
end

function [m,ft] = findfaceTension( m )
    [ft,faceNormals] = faceTension( m );
%     [vf,lengths] = invertIndexArray( m.FEconnectivity.faces, getNumberOfVertexes(m), 'cell' );
%     vt = zeros( getNumberOfVertexes(m), 1 );
%     for i=1:length(vf)
%         vt(i) = mean(ft(vf{i}));
%     end
%     v_tension_i = FindMorphogenIndex2( m, 'v_tension' );
%     if v_tension_i ~= 0
%         m.morphogens(:,v_tension_i) = vt;
%     end

    [~,v_level_p] = getMgenLevels( m, 'V_LEVEL' );
    v_tension_i = FindMorphogenIndex2( m, 'v_tension' );
    
    facelevels = v_level_p( m.FEconnectivity.faces );
    horizFaces = all( facelevels==facelevels(:,1), 2 );
    horizFaceVxs = m.FEconnectivity.faces(horizFaces,:);
    vt = maxArray( horizFaceVxs, repmat( ft(horizFaces), 1, 3 ), getNumberOfVertexes(m) );
    m.morphogens(:,v_tension_i) = vt;
end

function [m,result] = GFtbox_Postiterate_Callback( m ) %#ok<DEFNU>
    result = [];
    
    if getModelOption( m, 'crackbytension' )
        [m,ft] = findfaceTension( m );
        ft = ft * expDecay( m.globalProps.timestep, getModelOption( m, 'strainretention' ) );
    
        breakingstress = getModelOption( m, 'breakingstress' );
        v_fracture_i = FindMorphogenIndex2( m, 'v_fracture' );
        v_tension_i = FindMorphogenIndex2( m, 'v_tension' );
        if (v_fracture_i ~= 0) && ~isempty(breakingstress)
            
            m.morphogens(:,v_fracture_i) = double( m.morphogens(:,v_tension_i) > getModelOption( m, 'breakingstress' ) );
        end
    
        [~,s_weakness_p] = getMgenLevels( m, 'S_WEAKNESS' );
        breakingFaces = ft > getModelOption( m, 'breakingstress' ) .* (1-perVertextoperFace( m, s_weakness_p ));
        timedFprintf( 1, 'Postiterate: %d candidate faces to split.\n', sum(breakingFaces) );
        [m,splitVxs] = makeCrack( m, [], breakingFaces );
        timedFprintf( 1, 'Postiterate: %d vertexes split.\n', length(splitVxs) );
        m.userdata.numsplits( Steps(m)+1 ) = length(splitVxs);
        id_split_i = FindMorphogenIndex2( m, 'id_split' );
        m.morphogens(splitVxs,id_split_i) = 1;
    end
    xxxx = 1;
end

function [m,result] = GFtbox_Preplot_Callback( m, theaxes ) %#ok<DEFNU>
    result = [];

    [m,ft] = findfaceTension( m );
    
%     breakingstress = getModelOption( m, 'breakingstress' );
%     v_fracture_i = FindMorphogenIndex2( m, 'v_fracture' );
%     if (v_fracture_i ~= 0) && ~isempty(breakingstress)
%         m.morphogens(:,v_fracture_i) = double( vt > getModelOption( m, 'breakingstress' ) );
%     end
end

function [m,result] = GFtbox_Postplot_Callback( m, theaxes ) %#ok<DEFNU>
    result = [];
    % Your code here.
end

% function crackvxs = findRandomCracks( m, n, r, margin, extension )
% % N: number of rings of vertexes.
% % R: number of cracks to make.
% % MARGIN: number of planes of vertexes to ignore at top and bottom.
% % EXTENSION: number of steps to enlarge the crack by.
% % Some other parameters are hard-wired here:
% % minCrackWidth and maxCrackWidth: the min and max angles of the cracks.
% % numsplitvxrings: the depth of the crack in terms of vertex rings.
% %
% % The result is a boolean map of which vertexes have been selected.
% 
%     crackRings = selectRings( n-2*margin, r ) + margin;
%     crackAngles = selectCrackDirections( length(crackRings) );
%     
%     ringheights = m.userdata.zmin + (crackRings-1)*(m.userdata.height/(m.userdata.numvxtiers-1));
%     numsplitvxrings = 1;
%     crackvxs = false( getNumberOfVertexes(m), 1 );
%     minCrackWidth = pi/4;
%     maxCrackWidth = pi/2;
%     crackWidths = minCrackWidth + (maxCrackWidth-minCrackWidth)*rand(size(crackAngles));
%     for i=1:length(crackRings)
%         thisVxRing = abs(m.FEnodes(:,3) - ringheights(i)) < 0.001;
%         ca = crackAngles(i);
%         angleok = angleBetween( ca-crackWidths(i)/2, m.userdata.angles, ca+crackWidths(i)/2 );
%         radiusok = m.userdata.radii >= m.userdata.radius * (1 - (numsplitvxrings-1)/m.userdata.numelementrings - 0.001);
%         crackvxs = crackvxs | (thisVxRing & angleok & radiusok);
%     end
%     
%     if extension > 0
%         [~,v_level_p] = getMgenLevels( m, 'V_LEVEL' );
%         v_level_ends = v_level_p( m.FEconnectivity.edgeends );
%         planeedges = (v_level_ends(:,1) > 0) ...
%             & (v_level_ends(:,1) == v_level_ends(:,2));
%         ee = m.FEconnectivity.edgeends( planeedges, : );
%         for i=1:extension
%             extedges = any( crackvxs( ee ), 2 );
%             crackvxs( ee(extedges,:) ) = true;
%         end
%     end
% end

function crackvxs = minimalRandomCracks( m, n, r, margin, extension )
% N: number of rings of vertexes.
% R: number of cracks to make.
% MARGIN: number of planes of vertexes to ignore at top and bottom.
% EXTENSION: number of steps to enlarge the crack by.

    [~,v_level_p] = getMgenLevels( m, 'V_LEVEL' );
    [~,v_ring_p] = getMgenLevels( m, 'V_RING' );
    crackRings = selectRings( n-2*margin, r ) + margin;
    crackAngles = selectCrackDirections( length(crackRings) );
    
    crackvxs = false( getNumberOfVertexes(m), 1 );
    for i=1:length(crackRings)
        eligibleVxs = find( (v_level_p==crackRings(i)) & (v_ring_p==m.userdata.numelementrings) );
        ringAngles = m.userdata.angles( eligibleVxs );
        ca = crackAngles(i);
        [ang,angi] = nearestAngleTo( ca, ringAngles );
        chosenVx = eligibleVxs(angi);
        crackvxs(chosenVx) = true;
    end
    
    if extension > 0
        v_level_ends = v_level_p( m.FEconnectivity.edgeends );
        planeedges = (v_level_ends(:,1) > 0) ...
            & (v_level_ends(:,1) == v_level_ends(:,2));
        ee = m.FEconnectivity.edgeends( planeedges, : );
        for i=1:extension
            extedges = any( crackvxs( ee ), 2 );
            crackvxs( ee(extedges,:) ) = true;
        end
    end
end

function [ang,angi] = nearestAngleTo( a, angles )
% Find the member of ANGLES which is closest to the angle A,
% in circular measure, and return the angle and its index.

    a1 = mod( angles-a+pi, 2*pi ) - pi;
    [~,angi] = min(abs(a1));
    ang = angles(angi);
end

function ab = angleBetween( mina, a, maxa )
    ab = mod( a-mina, 2*pi ) < mod( maxa-mina, 2*pi );
end

function sr = selectRings( n, r )
    a = [ false randsubset( r, n-r ) ];
    b = ones(1,r+1);
    b(a) = b(a)+1;
    c = cumsum(b);
    c(end) = [];
    if (c(end)<n) && (rand(1) < 0.5)
        c = c+1;
    end
    sr = c;
end

function ca = selectCrackDirections( n )
    ca = mod( cumsum( rand(n,1) * pi + pi/2 ), 2*pi );
end

function ff = smoothField( m, numiters, range )
    ff = rand( getNumberOfVertexes(m), 1 );
    for i=1:numiters
        [ff1,nn] = averageArray( m.FEconnectivity.edgeends, ff(m.FEconnectivity.edgeends(:,[2 1])), size(ff) );
        ff = ff1;
    end
    ff = ff-min(ff);
    ff = ff * (range(2)-range(1))/max(ff);
    ff = ff + range(1);
end

function m = setProjectOptions( m )
    if m.globalProps.IFsetsoptions
        m = setUpModelOptions( m, ...
            'corerings', [], 5, ...
            'coreradius', [], 5/6 - 0.0004, ...
            'innergrowthrate', [], 0.04, ...
            'outergrowthrate', [], 0, ...  Set this to 0.04 for wild type model, and zero for dwarf model
            'innerhorizgrowthrate', [], 0, ...
            'outerhorizgrowthrate', [], 0, ...
            'fixoutervertexes', { 'z', 'xyz', '' }, '', ...
            'outerbulkmodulus', {}, 4, ... 
            'crackbytension', {}, true, ...
            'crackheight', {}, 0, ...
            'crackradius', {}, 2.5*4/6-0.0001, ...
            'crackstartangle', {}, -2*pi/3 - 0.0001, ...
            'crackendangle', {}, -pi/3 + 0.0001, ...
            'makeinitcracks', {}, false, ...
            'weakness', [], 0.4, ...
            'diffusiontime', [], 2, ...
            'breakingstress', {}, 0.07, ... 
            'strainretention', {}, 0.5 ...
        );
    end
end


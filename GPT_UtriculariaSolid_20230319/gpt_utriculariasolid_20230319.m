function [m,result] = gpt_utriculariasolid_20230319( m, varargin )
%[m,result] = gpt_utriculariasolid_20230319( m, varargin )
%   Morphogen interaction function.
%   Written at 2023-03-19 16:00:02.
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
        m = setProjectOptions( m );
        m = leaf_plotoptions( m, 'gradientoffset', 0.5, 'edgesharpness', pi/3 );
        m.outputcolors.residualgrowth = [ 0 0 1; 1 0 0 ]; % Red for tension, blue for compression.
        m = leaf_plotpriority( m, {'id_epidermis','id_blades','id_axis'}, [2 1 2], 0.5, 'type', 'morphogen' );
    end
    OPTIONS = getModelOptions( m );
    printModelOptions( m );
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
    [id_fixed_i,id_fixed_p,id_fixed_a,id_fixed_l] = getMgenLevels( m, 'ID_FIXED' );  %#ok<ASGLU>
    [id_axis_i,id_axis_p,id_axis_a,id_axis_l] = getMgenLevels( m, 'ID_AXIS' );  %#ok<ASGLU>
    [id_blades_i,id_blades_p,id_blades_a,id_blades_l] = getMgenLevels( m, 'ID_BLADES' );  %#ok<ASGLU>
    [id_epidermis_i,id_epidermis_p,id_epidermis_a,id_epidermis_l] = getMgenLevels( m, 'ID_EPIDERMIS' );  %#ok<ASGLU>
    [id_growth_i,id_growth_p,id_growth_a,id_growth_l] = getMgenLevels( m, 'ID_GROWTH' );  %#ok<ASGLU>
    [id_strainret_i,id_strainret_p,id_strainret_a,id_strainret_l] = getMgenLevels( m, 'ID_STRAINRET' );  %#ok<ASGLU>
    [v_stiffness_i,v_stiffness_p,v_stiffness_a,v_stiffness_l] = getMgenLevels( m, 'V_STIFFNESS' );  %#ok<ASGLU>
    [id_top_i,id_top_p,id_top_a,id_top_l] = getMgenLevels( m, 'ID_TOP' );  %#ok<ASGLU>
    [id_bottom_i,id_bottom_p,id_bottom_a,id_bottom_l] = getMgenLevels( m, 'ID_BOTTOM' );  %#ok<ASGLU>
    [id_endmargin_i,id_endmargin_p,id_endmargin_a,id_endmargin_l] = getMgenLevels( m, 'ID_ENDMARGIN' );  %#ok<ASGLU>
    [v_outside_i,v_outside_p,v_outside_a,v_outside_l] = getMgenLevels( m, 'V_OUTSIDE' );  %#ok<ASGLU>

% Mesh type: volumetric
%          FEtype: 84.00000081
%        axisdivs: 20
%      circumdivs: 18
%       coneangle: 0
%         dealign: 0
%          hollow: 0
%       innerdivs: 0
%             new: 1
%        position: 0
%           rings: 3
%            size: 1.00000012

%            Morphogen    Diffusion   Decay   Dilution   Mutant
%            --------------------------------------------------
%                 KPAR         ----    ----       ----     ----
%                KPAR2         ----    ----       ----     ----
%                 KPER         ----    ----       ----     ----
%                  POL         ----    ----       ----     ----
%                 POL2         ----    ----       ----     ----
%             ID_FIXED         ----    ----       ----     ----
%              ID_AXIS         ----    ----       ----     ----
%            ID_BLADES         ----    ----       ----     ----
%         ID_EPIDERMIS         ----    ----       ----     ----
%            ID_GROWTH         ----    ----       ----     ----
%         ID_STRAINRET         ----    ----       ----     ----
%          V_STIFFNESS         ----    ----       ----     ----
%               ID_TOP         ----    ----       ----     ----
%            ID_BOTTOM         ----    ----       ----     ----
%         ID_ENDMARGIN         ----    ----       ----     ----
%            V_OUTSIDE         ----    ----       ----     ----


%%% USER CODE: MORPHOGEN INTERACTIONS

% In this section you may modify the mesh in any way that does not
% alter the set of nodes.

    if (Steps(m)==0) && m.globalDynamicProps.doinit
        % Put any code here that should only be performed at the start of
        % the simulation.
        
        m.sharpedges = findSharpEdges( m, 1.22 );
        
        m = leaf_setmorphogenrole( m, 'id_strainret', 'strainret' );
        id_strainret_p=0.5; ...set strainretention to 0.5

        % These constants must agree with the dimensions of the mesh.
        AXIS_RADIUS = 0.5;
        EPIDERMIS_THICKNESS = AXIS_RADIUS*2/3;
        
        radii = sqrt( sum( m.FEnodes(:,[1 2]).^2, 2 ) );
        min_z = min( m.FEnodes(:,3) );
        max_z = max( m.FEnodes(:,3) );
                
        ABSTOL = EPIDERMIS_THICKNESS * 0.1;
        RELTOL = 0.001;
        id_top_p = m.FEnodes(:,3) >= max_z - ABSTOL;
        id_bottom_p = m.FEnodes(:,3) <= min_z + ABSTOL;
        axisnodes = radii <= AXIS_RADIUS + ABSTOL;
        id_axis_p = double( axisnodes );
        epidermisnodes = radii >= max(radii) - EPIDERMIS_THICKNESS - ABSTOL;
        id_epidermis_p = double( epidermisnodes );
        id_blades_p = 1 - max( id_axis_p, id_epidermis_p );
        v_outside_p = double( radii >= max(radii) - EPIDERMIS_THICKNESS + ABSTOL );
        
        % Don't do the following for the solid model.
        % Find the edges that join an epidermal vertex to a blades vertex.
%         epiEdges = any( id_epidermis_p( m.FEconnectivity.edgeends )==1, 2 ); % Edges with at least one end in the epidermis.
%         bladeEdges = any( id_blades_p( m.FEconnectivity.edgeends )==1, 2 ); % Edges with at least one end in the blades.
%         bladeEpiBoundaryEdges = epiEdges & bladeEdges; % Edges with one epidermal end and one blade end.
        % The vertexes where the blades join the epidermis are the epidermal ends of those edges.
%         bladeEpiBoundaryVxs = unique( m.FEconnectivity.edgeends( bladeEpiBoundaryEdges, : ) );
%         bladeEpiBoundaryVxs( id_blades_p(bladeEpiBoundaryVxs)==1 ) = [];
%         % Transfer these vertexes from id_blades to id_epidermis.
%         id_epidermis_p( bladeEpiBoundaryVxs ) = 0;
%         id_blades_p( bladeEpiBoundaryVxs ) = 1;
        
        if ~isempty( OPTIONS.fixoutervertexes )
            discRadius = max(radii);
            outerVxs = radii >= discRadius * (1 - RELTOL);
            switch OPTIONS.fixoutervertexes
                case 'z'
                    m.fixedDFmap( outerVxs, 3 ) = true;
                case 'xyz'
                    m.fixedDFmap( outerVxs, : ) = true;
            end
        end
        if ~isempty( OPTIONS.fixendvertexes )
            endVxs = abs( m.FEnodes(:,3) ) > max_z - ABSTOL;
            switch OPTIONS.fixendvertexes
                case 'z'
                    m.fixedDFmap( endVxs, 3 ) = true;
                case 'xyz'
                    m.fixedDFmap( endVxs, : ) = true;
            end
        end
        id_fixed_p = double( any( m.fixedDFmap, 2 ) );
        
        switch OPTIONS.growthregion
            case 'core'
                id_growth_p = id_axis_p;
            case 'dwarf'
                id_growth_p = max( id_axis_p, id_blades_p );
            case 'wt'
                id_growth_p(:) = 1;
            case 'blades'
                id_growth_p = id_blades_p;
            otherwise
                id_growth_p(:) = 0;
        end
        
        if OPTIONS.flatends
            % Constrain the top and bottom surfaces to remain flat, but
            % otherwise free to move.
            topZdfs = 3 * find( id_top_p==1 ); % Multiplying by 3 to get the Z d.o.f.
            bottomZdfs = 3 * find( id_bottom_p==1 ); % Multiplying by 3 to get the Z d.o.f.
            m.globalDynamicProps.stitchDFsets = { topZdfs, bottomZdfs };
        else
            m.globalDynamicProps.stitchDFsets = {};
        end
        
        if OPTIONS.endmargin > 0
            endmarginVxs = max(m.FEnodes(:,3)) - abs( m.FEnodes(:,3) ) < OPTIONS.endmargin;
            id_growth_p( endmarginVxs ) = 0;
            id_endmargin_p = double( endmarginVxs );
        else
            id_endmargin_p(:) = 0;
        end
        
        % Set epidermal stiffness.
        m.cellbulkmodulus(:) = 1;
        if OPTIONS.epidermalstiffness ~= 1
            epidermalElements = perVertextoperFE( m, id_epidermis_p, 'max' ) == 1;
            m.cellbulkmodulus( epidermalElements ) = OPTIONS.epidermalstiffness;
        end
        if OPTIONS.endmarginstiffness ~= 1
            endmarginElements = perVertextoperFE( m, id_endmargin_p, 'max' ) == 1;
            m.cellbulkmodulus( endmarginElements ) = OPTIONS.endmarginstiffness;
        end
        v_stiffness_p = perFEtoperVertex( m, m.cellbulkmodulus, 'max' );

        % Vertical gradient of pol.
        pol_p = m.FEnodes(:,3);
        pol_p = pol_p - min(pol_p);
        pol_p = pol_p/max(pol_p);
    
        % Radial gradient of pol2
        pol2_p = radii;
    end
    
    if Steps(m)==1
        % Apply a smooth random perturbation to the mesh to break symmetry.
        RANDOMISE = true;
        if RANDOMISE
            radii = sqrt( sum( m.FEnodes(:,[1 2]).^2, 2 ) );
            [zz,p] = sort( m.FEnodes(:,3) );
            [clumpindex,clumpwidths,~] = clumpValues( zz, 0.001 );
            randXY0 = randn( length(clumpwidths), 2 );
            f = exp( -(linspace( -2, 2, 21 ).^2 )/2 )';
            randXY(:,1) = imfilter( randXY0(:,1), f, 'circular', 'conv' );
            randXY(:,2) = imfilter( randXY0(:,2), f, 'circular', 'conv' );
            randXY = randXY / std(randXY(:));
            randXY(randXY > 3) = 3;
            randXY(randXY <= -3) = -3;
            randXY = randXY * (max(radii) * 0.01);

            offset = [ randXY( clumpindex, : ), zeros( getNumberOfVertexes(m), 1 ) ];
            m.FEnodes(p,:) = m.FEnodes(p,:) + offset;
        end
        
        switch getModelOption( m, 'variant' )
            case 'AXIAL'
                axialgrowthrate = OPTIONS.growthrate;
                radialgrowthrate = 0;
            case 'RADIAL'
                axialgrowthrate = 0;
                radialgrowthrate = OPTIONS.growthrate;
            otherwise
                axialgrowthrate = 0;
                radialgrowthrate = 0;
        end
        kpar_p = axialgrowthrate * id_growth_p;
        kpar2_p(:) = radialgrowthrate * id_blades_p;
        kper_p(:) = 0;
    end
%%% END OF USER CODE: MORPHOGEN INTERACTIONS

%%% SECTION 3: INSTALLING MODIFIED VALUES BACK INTO MESH STRUCTURE
%%% AUTOMATICALLY GENERATED CODE: DO NOT EDIT.

    m.morphogens(:,kpar_i) = kpar_p;
    m.morphogens(:,kpar2_i) = kpar2_p;
    m.morphogens(:,kper_i) = kper_p;
    m.morphogens(:,pol_i) = pol_p;
    m.morphogens(:,pol2_i) = pol2_p;
    m.morphogens(:,id_fixed_i) = id_fixed_p;
    m.morphogens(:,id_axis_i) = id_axis_p;
    m.morphogens(:,id_blades_i) = id_blades_p;
    m.morphogens(:,id_epidermis_i) = id_epidermis_p;
    m.morphogens(:,id_growth_i) = id_growth_p;
    m.morphogens(:,id_strainret_i) = id_strainret_p;
    m.morphogens(:,v_stiffness_i) = v_stiffness_p;
    m.morphogens(:,id_top_i) = id_top_p;
    m.morphogens(:,id_bottom_i) = id_bottom_p;
    m.morphogens(:,id_endmargin_i) = id_endmargin_p;
    m.morphogens(:,v_outside_i) = v_outside_p;

%%% USER CODE: FINALISATION

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

function [m,result] = GFtbox_Precelldivision_Callback( m, ci ) %#ok<INUSD>
    result = [];
    % Your code here.

    % If a nonempty result is to be returned, it should be a struct
    % with fields result.divide, result.dividepoint, and result.perpendicular.
end

function [m,result] = GFtbox_Postcelldivision_Callback( m, ci, cei, newci, newcei, oe1, oe2, ne1, ne2, ne3 ) %#ok<INUSD>
    result = [];
    % Your code here.
end

function [m,result] = GFtbox_Postiterate_Callback( m )
    result = [];
    % Your code here.
end

function [m,result] = GFtbox_Preplot_Callback( m, theaxes ) %#ok<INUSD>
    result = [];
    % Your code here.
end

function [m,result] = GFtbox_Postplot_Callback( m, theaxes ) %#ok<INUSD>
    result = [];
    % Your code here.
end

function m = setProjectOptions( m )
    if m.globalProps.IFsetsoptions
        m = setUpModelOptions( m, ...
            'variant', { 'AXIAL', 'RADIAL' }, 'AXIAL', ...
            'growthrate', [], 0.04, ...
            'growthregion', { 'core', 'dwarf', 'wt', 'blades' }, 'dwarf', ...
            ... % 'dwarf': Utric Dwarf. Growth in the core and the blades.
            ... % 'wt': wild type. Uniform growth everywhere.
            ... % 'core': Growth only in the central cylinder.
            ... % 'blades' Growth only in the blades.
            ... % The AXIAL variant can be used with any growth region.
            ... % The RADIAL variant is intended to be used only with the 'blades' growth region.
            ...
            'epidermalstiffness', {}, 2, ... % 1 means the same as everywhere else.
            'endmarginstiffness', {}, 1, ... % 1 means the same as everywhere else.
            'endmargin', [], 0, ... % No growth within this distance of the initial top or bottom. Default should be 1 if flatends false, 0 if true.
            'flatends', [], true, ... % Constrain the top and bottom faces to remain exactly flat in the XY plane.
            'fixoutervertexes', { 'z', 'xyz', '' }, '', ... % Fix degrees of freedom of the vertexes on the outer cylindrical surface.
            'fixendvertexes', { 'z', 'xyz', '' }, '', ... % Fix degrees of freedom of the vertexes on the top and bottom surfaces.
            ...
            'cluster_expt_id', {}, '', ... % Used for storing information about cluster runs.
            'cluster_expt_variant', {}, '', ... % Used for storing information about cluster runs.
            'cluster_expt_repetition', {}, '', ... % Used for storing information about cluster runs.
            'comment', {}, '' ... % Used for storing information about cluster runs.
        );
    end
end



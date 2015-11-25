function [u, dirs, mesh] = platonicSolid(shape)
%PLATONICSOLID Generates platonic solid vertices and faces
% PLATONICSOLID Generates the five platonic solids, centered at the origin,
% with their vertices being at unit radius. The cartesian coordinates, the
% directions (in [azi elev] convention), or the mesh structure with
% vertices and faces can be returned. The solid to be generated is
% specified by the string 'shape'.
%
% Inputs:
%   shape:  'tetra' or 'tetrahedron'
%           'hexa' or 'hexahedron' or 'cube'
%           'octa' or 'octahedron'
%           'dodeca' or 'dodecahedron'
%           'icosa' or 'icosahedron'
%           string to get the respective shape, capital or lowercase
%
% Outputs:
%   u:      matrix of vertex coordinates
%   dirs:   matrix of spherical coordinates in [azi elev] style
%   mesh:   structure combining vertices and faces, for plotting 3D plots,
%           or generating dense meshes by subdividing faces and
%           re-triangulating
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Archontis Politis, 15/11/2015
%   archontis.politis@aalto.fi
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

switch lower(shape)
    case {'tetra','tetrahedron'}
        u = [1     1     1
            1    -1    -1
            -1     1    -1
            -1    -1     1 ]/sqrt(3);
        
        f = [1     2     3
            1     3     4
            1     4     2
            2     4     3 ];
        
    case {'cube','hexa','hexahedron'}
        u = [1     1     1
            1     1    -1
            1    -1     1
            1    -1    -1
            -1     1     1
            -1     1    -1
            -1    -1     1
            -1    -1    -1 ]/sqrt(3);
        
        f = [2     4     3     1
            8     6     5     7
            6     2     1     5
            4     8     7     3
            3     7     5     1
            6     8     4     2];
        
    case {'octa','octahedron'}
        u = [1     0     0
            -1     0     0
            0     1     0
            0    -1     0
            0     0     1
            0     0    -1 ];
        
        f = [1     3     5
            1     4     6
            1     5     4
            1     6     3
            2     3     6
            2     4     5
            2     5     3
            2     6     4];
        
    case {'dodeca','dodecahedron'}
        gr = (1 + sqrt(5))/2; % golden ratio
        
        u = [1     1      -1
            -1     1      -1
            -1     1       1
            -1    -1       1
            1    -1       1
            1    -1      -1
            -1    -1      -1
            1     1       1
            0     1/gr   -gr
            0    -1/gr    gr
            0    -1/gr   -gr
            0     1/gr    gr
            1/gr  gr      0
            -1/gr  gr      0
            1/gr -gr      0
            -1/gr -gr      0
            -gr    0      -1/gr
            -gr    0       1/gr
            gr    0       1/gr
            gr    0      -1/gr ] / sqrt(3);
        
        f = [2     3       5       4       1
            17    14      15      2       1
            18    16      11      14      17
            1     4       6       18      17
            15    13      19      3       2
            14    11      12      13      15
            13    12      10      20      19
            3     19      20      7       5
            4     5       7       8       6
            10    9       8       7       20
            18    6       8       9       16
            16    9       10      12      11 ];
                
    case {'icosa','icosahedron'}
        gr  = (1 + sqrt(5)) / 2; % golden ratio
        
        u = [0    gr    1
            0   -gr    1
            0    gr   -1
            0   -gr   -1
            1    0    gr
            1    0   -gr
            -1    0    gr
            -1    0   -gr
            gr    1    0
            gr   -1    0
            -gr    1    0
            -gr   -1    0 ] / sqrt(1+gr^2);
        
        f = [1     3    11
            1     5     9
            1     7     5
            1     9     3
            1    11     7
            2     4    10
            2     5     7
            2     7    12
            2    10     5
            2    12     4
            3     6     8
            3     8    11
            3     9     6
            4     6    10
            4     8     6
            4    12     8
            5    10     9
            6     9    10
            7    11    12
            8    12    11 ];
        
end

% return directions and mesh with vertices/faces
[dirs(:,1), dirs(:,2)] = cart2sph(u(:,1), u(:,2), u(:,3));
mesh.vertices = u;
mesh.faces = f;

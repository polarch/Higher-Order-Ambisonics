function [u, dirs, mesh] = platonicSolid(shape)
%PLATONICSOLID Generates platonic solid vertices and faces
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Archontis Politis, 15/11/2015
%   archontis.politis@aalto.fi
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

switch lower(shape)
    case {'tetra','tetrahedron'}
        u = [ 1     1     1
              1    -1    -1
             -1     1    -1
             -1    -1     1 ]/sqrt(3);
        f = [ 1    2    3
              1    4    2
              3    2    4
              4    1    3 ];
    case {'cube','hexa','hexahedron'}
        u = [ 1     1     1
              1     1    -1
              1    -1     1
              1    -1    -1
             -1     1     1
             -1     1    -1
             -1    -1     1
             -1    -1    -1 ]/sqrt(3);        
        f =  [ 2    4    3    1
               8    6    5    7
               6    2    1    5
               4    8    7    3
               3    7    5    1
               6    8    4    2];
    case {'octa','octahedron'}
        u = [ 1     0     0
               -1     0     0
                0     1     0
                0    -1     0
                0     0     1
                0     0    -1 ];
        f = [ 1    5    3
                3    5    2
                2    5    4
                4    5    1
                1    3    6
                3    2    6
                2    4    6
                4    1    6 ];
    case {'dodec','dodecahedron'}
        % Golden mean
        gm = (sqrt(5)+1)/2;
        
        u = [ 0     1/gm   -gm
              1     1      -1
              1/gm  gm      0
             -1     1      -1
             -1/gm  gm      0
             -gm    0      -1/gm
             -1     1       1
             -gm    0       1/gm
             -1    -1       1
              0    -1/gm    gm
              1/gm -gm      0
              1    -1       1
              gm    0       1/gm
              1    -1      -1
              gm    0      -1/gm
             -1/gm -gm      0
              0    -1/gm   -gm
             -1    -1      -1
              1     1       1
              0     1/gm    gm ] / sqrt(3);
        
        f =  [ 2     3     5     4     1
                   17    14    15     2     1
                   18    16    11    14    17
                   1     4     6    18    17
                   15    13    19     3     2
                   14    11    12    13    15
                   13    12    10    20    19
                   3    19    20     7     5
                   4     5     7     8     6
                   10     9     8     7    20
                   18     6     8     9    16
                   16     9    10    12    11 ];
        
    case {'icosa','icosahedron'}
        
        % Golden mean
        gm  = (sqrt(5)+1) / 2;
        tau = gm / sqrt(1 + gm^2);
        one = 1 / sqrt(1 + gm^2);
        
        u = [  tau     one       0
              -tau     one       0
              -tau    -one       0
               tau    -one       0
               one       0     tau
               one       0    -tau
              -one       0    -tau
              -one       0     tau
               0     tau     one
               0    -tau     one
               0    -tau    -one
               0     tau    -one ];
        
        f = [  5     8     9
            5    10     8
            6    12     7
            6     7    11
            1     4     5
            1     6     4
            3     2     8
            3     7     2
            9    12     1
            9     2    12
            10     4    11
            10    11     3
            9     1     5
            12     6     1
            5     4    10
            6    11     4
            8     2     9
            7    12     2
            8    10     3
            7     3    11 ];
        
end

[dirs(:,1), dirs(:,2)] = cart2sph(u(:,1), u(:,2), u(:,3));
mesh.vertices = u;
mesh.faces = f;

function C = cr_coarsening(ia,ja,aa,nu,alpha)
% CR_COARSENING Compatible relaxation based coarsening
%   C = CR_COARSENING(IA,JA,AA,NU,ALPHA)
%
%   See also amg_coarsening, cdata.

%   Written by Minho Park
%   Department of Applied Mathematics
%   University of Colorado at Boulder
%   parkmh@colorado.edu
%   $Date: 2009/02/12 $
%   Update : 

DEBUG = 0;

if DEBUG
    fid = fopen('Coarsening.dbg','w');
    fprintf(fid,' *** First coloring pass ***\n');
end

m = length(ia)-1;   % Size of matrix A

% Memory allocation
C = zeros(1,m);     % Coarse grid point
U = zeros(1,m);     % Coarse grid candidate
w = zeros(1,m);     % weight
second_coloring_pt = zeros(1,m);

uindex = 0;
cg_ratio = 0;   % coarsening ratio
second_index = 1;

% Random initial
x = 1 + rand(m,1)/4;


% From the first to one sweep before the last
gs(ia,ja,aa,x,zeros(m,1),1,nu-1);
bottom = sqrt(x'*mm(ia,ja,aa,x));

% The last iteration sweep
gs(ia,ja,aa,x,zeros(m,1),1,1);
top = sqrt(x'*mm(ia,ja,aa,x));

% stoppint criterion
rho_cr = top/bottom;    

while rho_cr >= alpha
    max_x = max(abs(x));
    
    % find coarse grid candidate U
    for i = 1:m
        if ~C(i)
            if abs(x(i)) >= (1 - rho_cr)*max_x
                U(i) = 1;
                uindex = uindex + 1;
            else
                second_coloring_pt(i) = 1;
                second_index = second_index + 1;
            end
        end
    end
    
    % calculate weight
    for i = 1:m
        if U(i)
            for j = ia(i) : ia(i+1)-1
                if U(ja(j))
                    w(i) = w(i) + 1;
                end
            end
        end
    end 
    
    windex = 2;
    [sorted_w index_from] = sort(w,'descend');
    [dummy_index index_to] = sort(index_from);
    
    while uindex > 0
        if windex == m
            break
        end
        unassigned = U;
        
        % find first nonzeros
        while ~w(windex-1)
            windex = windex + 1;
        end
        
        C( index_from( windex - 1 ) ) = 1;
        U( index_from( windex - 1 ) ) = 0;
        uindex = uindex - 1;
        w( windex - 1 ) = 0;
        
        for j = ia(index_from(windex - 1)):ia(index_from(windex - 1) + 1) - 1
            if U(ja(j))
                unassigned(ja(j)) = 0;
            end
        end
        
        for j = ia(index_from(windex - 1)):ia(index_from(windex - 1) + 1) - 1
            if U(ja(j))
                U(ja(j)) = 0;
                uindex = uindex - 1;
                w( index_to( ja(j) ) ) = 0;
                for k = ia(ja(j)) : ia(ja(j)+ 1) - 1
                    if U(ja(k)) && unassigned(ja(k))
                        w(index_to(ja(k))) = w(index_to(ja(k))) + 1;
                        for l = ia(ja(k)):ia(ja(k)+ 1)-1
                            if U(ja(l))&&unassigned(ja(l))
                                w(index_to(ja(l))) = w(index_to(ja(l))) + 1;
                                if w(index_to(ja(l))) < 0
                                    w(index_to(ja(l))) = 0;
                                end
                            end
                        end
                    end
                end
            end
        end
    end
    
    cg_ratio = sum(C)/m;
    x = 1 + rand(m,1)/4;  
    x(C == 1) = 0;
    
    % From the first to one sweep before the last
    cr(ia,ja,aa,x,zeros(m,1),1.0,nu - 1,C);
    bottom = sqrt(x'*mm(ia,ja,aa,x));

    % The last iteration sweep
    cr(ia,ja,aa,x,zeros(m,1),1.0,1,C);
    top = sqrt(x'*mm(ia,ja,aa,x));
    
    if bottom == 0
        fprintf('it breaks here \n')        
        break
    end
    
    rho_cr = top/bottom;
    
end


% Second coloring
int connection2C;
if second_index > 0
    for i = 1:m
        if second_coloring_pt(i)
            connection2C = 0;
            if ia(i)<ia(i+1)-1
            for j = ia(i) : ia(i+1)-1
                if C(ja(j))
                    connection2C = 1;
                end
            end
            end
            if connection2C == 0
                C(i) = 1;
            end
        end
    end
end

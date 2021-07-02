% v.0.1.2 - New amgbench.m
% TODO test Wcycle result
% TODO test complex valued matrix
% TODO test adaptive AMG and rBAMG

clear all

today = date;
logfile = [datestr(datenum(today),29) '.bench'];
amgver = amgversion;
fid = fopen(logfile,'w');
texfid = fopen('amgbench.tex','w');
title = 'Report on AMG Performance';
author = 'Minho Park';
opt = amgset;
    
N = [64 128 256 512 1024];
maxTest = 1;
fprintf(fid,'***********************************\n');
fprintf(fid,'%20s\n','AMG Reviews');
fprintf(fid,'%25s\n',['Date : ' today]);
fprintf(fid,'%21s\n',['ver : ' amgver]);
fprintf(fid,'***********************************\n\n');

fprintf(texfid,'\\documentclass[10pt]{report}\n');
fprintf(texfid,'\\title{%s}\n',title);
fprintf(texfid,'\\author{%s}\n',author);
fprintf(texfid,'\\begin{document}\n');
fprintf(texfid,'\\maketitle\n');
fprintf(texfid,'\\tableofcontents\n');
for chap = 1 : 3
    fprintf(fid,'\n\n######################################################################\n');
    switch chap
        case 1
            fprintf(fid,'%d. AMG on 2D Poisson (5 point stencil)\n',chap);
            fprintf(fid,'--------------------------------------\n');
            stencil = [4 -1 -1 -1 -1];
            fprintf(texfid,'\\chapter{2D Poisson (5 point stencil)}\n');
        case 2
            fprintf(fid,'%d. AMG on 2D Poisson (9 point stencil)\n',chap);
            fprintf(fid,'--------------------------------------\n');
            stencil = [8 -1 -1 -1 -1 -1 -1 -1 -1];
            fprintf(texfid,'\\newpage\n');
            fprintf(texfid,'\\chapter{2D Poisson (9 point stencil)}\n');
        case 3
            fprintf(fid,'%d. AMG on anisotropic 2D Poisson (5 point stencil)\n',chap);
            fprintf(fid,'--------------------------------------\n');
            epsilon = .1;
            stencil = [2*(1+epsilon) -epsilon -epsilon -1 -1];
            fprintf(texfid,'\\newpage\n');
            fprintf(texfid,'\\chapter{Anisotropic 2D Poisson (9 point stencil)}\n');
    end
    for section = 1 : 2
        if section == 1
            fprintf(fid,'\n%d.%d. Classical AMG\n',chap,section);
            fprintf(fid,'------------------\n');
            opt = amgset(opt,'PreCond','off');
            fprintf(texfid,'\\section{Performance of AMG}\n');
        elseif section == 2
            fprintf(fid,'\n%d.%d. PCG with AMG\n',chap,section);
            fprintf(fid,'------------------\n');
            opt = amgset(opt,'PreCond','pcg');
            fprintf(texfid,'\\section{Performance of PCG}\n');
        end
        
        for subsection = 1 : maxTest
            fprintf(fid,'\n%d.%d.%d. N = %d\n',chap,section,subsection,N(subsection));
            fprintf(fid,'------------------\n');
            fprintf(texfid,'\\subsection{N = %d}\n',N(subsection));
            A = stencil2mat(stencil,N(subsection));
            opt = amgset(opt,'Coarsest',5+subsection);
            b = zeros(length(A),1);
            v = rand(size(b));
            profile clear;
            profile on;            
            [v amgout amgtime] = amg(A,v,b,opt);
            profile off;
            p = profile('info');
            profsave(p,[num2str(chap) '.' num2str(section) '.' num2str(subsection)]);
            allnnzs = 0; alldofs = 0;
            
            for subsubsection = 1 : 5
                switch subsubsection
                    case 1
                        fprintf(fid,'\n%d.%d.%d.%d. Operator Complexity\n',...
                            chap,section,subsection,subsubsection);
                        fprintf(fid,'-----------------------------\n');
                        fprintf(fid,'=================\n');
                        fprintf(fid,'LEVEL ||     NNZs\n');
                        fprintf(fid,'=================\n');
                        for i = 1 : amgout.cst
                            fprintf(fid,'%5d || %8d\n',i,amgout.nnz(i));
                            allnnzs = allnnzs + amgout.nnz(i);
                        end
                        fprintf(fid,'*** Operator Complexity : %2.2f [(nnz(1)+...+nnz(%d))/nnz(1)]\n',...
                            allnnzs/amgout.nnz(1),amgout.cst);                        
                    case 2
                        fprintf(fid,'\n%d.%d.%d.%d. Grid Complexity\n',...
                            chap,section,subsection,subsubsection);
                        fprintf(fid,'-------------------------\n');
                        fprintf(fid,'=================\n');
                        fprintf(fid,'LEVEL ||     DOFs\n');
                        fprintf(fid,'=================\n');
                        for i = 1 : amgout.cst
                            fprintf(fid,'%5d || %8d\n',i,amgout.nofdof(i));
                            alldofs = alldofs + amgout.nofdof(i);
                        end
                        fprintf(fid,'*** Grid Complexity     : %2.2f [(dof(1)+...+dof(%d))/dof(1)]\n',alldofs/amgout.nofdof(1),amgout.cst);                        
                    case 3
                        fprintf(fid,'\n%d.%d.%d.%d. Density (nnz per each row)\n',...
                            chap,section,subsection,subsubsection);
                        fprintf(fid,'------------------------------------\n');
                        fprintf(fid,'=================\n');
                        fprintf(fid,'LEVEL ||  Density\n');
                        fprintf(fid,'=================\n');
                        for i = 1 : amgout.cst
                            fprintf(fid,'%5d || %6.1f %% \n',i,amgout.nnz(i)/amgout.nofdof(i)^2*100);
                        end
                    case 4
                        fprintf(fid,'\n%d.%d.%d.%d. Convergence History\n',...
                            chap,section,subsection,subsubsection);
                        fprintf(fid,'-----------------------------\n');
                        fprintf(fid,'=====================================================\n');
                        fprintf(fid,'%5s ||%10s%6s%6s |%9s%6s%6s\n','CYCLE','RESIDUAL', 'CF', 'AVG','ERROR','CF','AVG');
                        fprintf(fid,'=====================================================\n');
                        fprintf(fid,'%3s   ||%10.2e%6s%6s |%9.2e%6s%6s\n','0',amgout.residual(1), '', '',amgout.errors(1),'','');
                        for i = 1 : amgout.amgcycle
                            fprintf(fid,'%3d   ||%10.2e%6.2f%6.2f |%9.2e%6.2f%6.2f\n',i,amgout.residual(i+1), amgout.residual(i+1)/amgout.residual(i),...
                                (amgout.residual(i+1)/amgout.residual(1))^(1/i),amgout.errors(i+1),amgout.errors(i+1)/amgout.errors(i),...
                                (amgout.errors(i+1)/amgout.errors(1))^(1/i));
                        end    
                    case 5
                        fprintf(fid,'\n%d.%d.%d.%d. Time\n',...
                            chap,section,subsection,subsubsection);
                        fprintf(fid,'----------------\n');
                        
                        fprintf(fid,'====================================================================\n');
                        fprintf(fid,'%5s || %14s %14s %14s %14s\n','LEVEL','S','C','P','RAP');
                        fprintf(fid,'====================================================================\n');
                        
                        setupTime = zeros(amgout.cst-1,1);
                        for i = 1 : amgout.cst - 1
                            setupTime(i) = amgtime.S(i)+ amgtime.C(i)+amgtime.P(i)+amgtime.RAP(i);
                            fprintf(fid,'%5d || %5.1e[%4.1f%%] %5.1e[%4.1f%%] %5.1e[%4.1f%%] %5.1e[%4.1f%%]\n',i,...
                                amgtime.S(i),amgtime.S(i)*100/setupTime(i),...
                                amgtime.C(i),amgtime.C(i)*100/setupTime(i),...
                                amgtime.P(i),amgtime.P(i)*100/setupTime(i),...
                                amgtime.RAP(i),amgtime.RAP(i)*100/setupTime(i));
                        end
                        fprintf(fid,'%5s || %5.1e[%4.1f%%] %5.1e[%4.1f%%] %5.1e[%4.1f%%] %5.1e[%4.1f%%]\n',...
                                'total',...
                                sum(amgtime.S),sum(amgtime.S)*100/sum(setupTime),...
                                sum(amgtime.C),sum(amgtime.C)*100/sum(setupTime),...
                                sum(amgtime.P),sum(amgtime.P)*100/sum(setupTime),...
                                sum(amgtime.RAP),sum(amgtime.RAP)*100/sum(setupTime));
                        fprintf(fid,'*** Total setup time : %.2e[%.2e](s)\n',sum(setupTime),amgtime.total);
                        fprintf(fid,'*** time for AMG cycles : %.2e(s) [avg : %.2e(s)]\n',amgtime.cycle,amgtime.cycle/amgout.amgcycle);
                        fprintf(fid,'*** AMG setup time is almost equivalent to %0.2f AMG cycles\n',sum(setupTime)*amgout.amgcycle/amgtime.cycle);
                    
                end
            end
            
            fprintf(texfid,'\\subsubsection{Complexties}\n');
            fprintf(texfid,'\\begin{table}[ht]\n');
            fprintf(texfid,'\\centering\n');
            fprintf(texfid,'\\begin{tabular}{|c||c|c|c|}\n');
            fprintf(texfid,'\\hline\n');
            fprintf(texfid,'Level & nnz & dof & density\\\\\\hline\n');
            fprintf(texfid,'\\hline\n');
            for i = 1 : amgout.cst
                fprintf(texfid,'%d & %d & %d & %4.2f\\\\\\hline\n',i,amgout.nnz(i),amgout.nofdof(i),amgout.nnz(i)/amgout.nofdof(i)^2*100);
            end
            fprintf(texfid,'\\hline\n');
            fprintf(texfid,' & Operator Complexity & Grid Complexity & \\\\\\hline\n');
            fprintf(texfid,' & %.2f & %.2f & \\\\\\hline\n',allnnzs/amgout.nnz(1),alldofs/amgout.nofdof(1));
            fprintf(texfid,'\\end{tabular}\n');
            fprintf(texfid,'\\end{table}\n');
            
            fprintf(texfid,'\\newpage\n');
            fprintf(texfid,'\\subsubsection{Convergence History}\n');
            fprintf(texfid,'\\begin{table}[ht]\n');
            fprintf(texfid,'\\centering\n');
            fprintf(texfid,'\\begin{tabular}{|c||c c c|c c c|}\n');
            fprintf(texfid,'\\hline\n');
            fprintf(texfid,'Cycle & Residual & cf & avg cf & error & cf & avg cf\\\\\\hline\n');
            fprintf(texfid,'\\hline\n');
            for i = 1 : amgout.amgcycle
                fprintf(texfid,'%d & %.1e & %.2f & %.2f & %.1e & %.2f & %.2f \\\\\\hline\n',...
                    i,amgout.residual(i+1), amgout.residual(i+1)/amgout.residual(i),...
                    (amgout.residual(i+1)/amgout.residual(1))^(1/i),amgout.errors(i+1),amgout.errors(i+1)/amgout.errors(i),...
                    (amgout.errors(i+1)/amgout.errors(1))^(1/i));
            end
            fprintf(texfid,'\\end{tabular}\n');
            fprintf(texfid,'\\end{table}\n');
            
            fprintf(texfid,'\\newpage\n');
            fprintf(texfid,'\\subsubsection{Time}\n');
            fprintf(texfid,'\\begin{table}[ht]\n');
            fprintf(texfid,'\\centering\n');
            fprintf(texfid,'\\begin{tabular}{|c||c c c c|}\n');
            fprintf(texfid,'\\hline\n');
            fprintf(texfid,'Level & S & C & P & RAP\\\\\\hline\n');
            fprintf(texfid,'\\hline\n');
            for i = 1 : amgout.cst - 1
                fprintf(texfid,'%d & %.1e[%4.1f\\%%] & %.1e[%4.1f\\%%] & %.1e[%4.1f\\%%] & %.1e[%4.1f\\%%] \\\\\\hline\n',i,...
                    amgtime.S(i),amgtime.S(i)*100/setupTime(i),...
                    amgtime.C(i),amgtime.C(i)*100/setupTime(i),...
                    amgtime.P(i),amgtime.P(i)*100/setupTime(i),...
                    amgtime.RAP(i),amgtime.RAP(i)*100/setupTime(i));
            end
            fprintf(texfid,'\\end{tabular}\n');
            fprintf(texfid,'\\end{table}\n');
        end
    end
end



fprintf(texfid,'\\end{document}\n');
fclose(texfid);
fclose(fid);

delete *.aux *.dvi *.log *.toc
fprintf('Creating amgbench.pdf\n')


system('latex -interaction=nonstopmode amgbench.tex','-echo');

system('pdflatex -interaction=nonstopmode amgbench.tex','-echo');
% system('xpdf amgbench.pdf');
clear


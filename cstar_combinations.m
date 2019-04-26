function [resultsX, qq] = cstar_combinations(cstar, step)
    % % All of the combinations of the properties used
    % cstar=[0.01 0.1 1.0 10.0];
    % dH=[20000 50000 80000 100000 150000 200000];
    % alp=[0.01 0.05 0.1 0.2 0.5 1.0];
    % % step taken for the mass fraction of each volatility bin tested
    % step=0.1;
    qq=0;
    % finds all combinations for volatility distribution if there are 4 bins
    if length(cstar)==4
        for ai=0:step:1
            for aj=0:step:(1.0-ai)
                for ak=0:step:(1.0-ai-aj)
                    qq=qq+1;
                    X_i(1)=ai;
                    X_i(2)=aj;
                    X_i(3)=ak;
                    X_i(4)=(1.0-ai-aj-ak);
                    X_a(qq,:)=X_i;
                    resultsX(qq, :) = X_i;
                end
            end
        end
        
        % finds all combinations for volatility distribution if there are 5 bins
    elseif length(cstar)==5
        for bi=0:step:1
            for bj=0:step:(1.0-bi)
                for bk=0:step:(1.0-bi-bj)
                    for bl=0:step:(1.0-bi-bj-bk)
                        qq=qq+1;
                        X_i(1)=bi;
                        X_i(2)=bj;
                        X_i(3)=bk;
                        X_i(4)=bl;
                        X_i(5)=(1.0-bi-bj-bk-bl);
                        X_a(qq,:)=X_i;
                        % %
                        resultsX(qq,:)= X_i;
                    end
                end
            end
        end
        % finds all combinations for volatility distribution if there are 6 bins
    elseif length(cstar)==6
        for ci=0:step:1
            for cj=0:step:(1.0-ci)
                for ck=0:step:(1.0-ci-cj)
                    for cl=0:step:(1.0-ci-cj-ck)
                        for cm=0:step:(1.0-ci-cj-ck-cl)
                            qq=qq+1;
                            X_i(1)=ci;
                            X_i(2)=cj;
                            X_i(3)=ck;
                            X_i(4)=cl;
                            X_i(5)=cm;
                            X_i(6)=(1.0-ci-cj-ck-cl-cm);
                            X_a(qq,:)=X_i;
                            % %
                            resultsX(qq,:)= X_i;
                        end
                    end
                end
            end
        end
        % finds all combinations for volatility distribution if there are 7 bins
    elseif length(cstar)==7
        for di=0:step:1
            for dj=0:step:(1.0-di)
                for dk=0:step:(1.0-di-dj)
                    for dl=0:step:(1.0-di-dj-dk)
                        for dm=0:step:(1.0-di-dj-dk-dl)
                            for dn=0:step:(1.0-di-dj-dk-dl-dm)
                                qq=qq+1;
                                X_i(1)=di;
                                X_i(2)=dj;
                                X_i(3)=dk;
                                X_i(4)=dl;
                                X_i(5)=dm;
                                X_i(6)=(dn);
                                X_i(7)=(1.0-di-dj-dk-dl-dm-dn);
                                X_a(qq,:)=X_i;
                                % %
                                resultsX(qq,:)= X_i;
                            end
                        end
                    end
                end
            end
        end
        % finds all combinations for volatility distribution if there are 8 bins
    elseif length(cstar)==8
        for ei=0:step:1
            for ej=0:step:(1.0-ei)
                for ek=0:step:(1.0-ei-ej)
                    for el=0:step:(1.0-ei-ej-ek)
                        for em=0:step:(1.0-ei-ej-ek-el)
                            for en=0:step:(1.0-ei-ej-ek-el-em)
                                for eo=0:step:(1.0-ei-ej-ek-el-em-en)
                                    qq=qq+1;
                                    X_i(1)=ei;
                                    X_i(2)=ej;
                                    X_i(3)=ek;
                                    X_i(4)=el;
                                    X_i(5)=em;
                                    X_i(6)=en;
                                    X_i(7)=eo;
                                    X_i(8)=(1.0-ei-ej-ek-el-em-en-eo);
                                    X_a(qq,:)=X_i;
                                    % %
                                    resultsX(qq,:)= X_i;
                                end
                            end
                        end
                    end
                end
            end
        end
    end
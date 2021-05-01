function [xf,ff,xVec,fVec] = newtonRaphsonND(f,df,x0,fTol,xTol,maxIt,figHandle)
%NEWTONRAPHSON1D Summary of this function goes here
%   Detailed explanation goes here
    
    [nrows,ncols] = size(x0);
    if nrows == 1 || ncols == 1
        %input is a vector, good
        if ncols == 1
            %input is a column vector, good
            n = nrows;
        elseif nrows == 1
            %input is a row vector, transpose it
            x0 = x0';
            n = ncols;
        end
    else
        error('x must be a vector, not a matrix.')
    end
        
    f0 = f(x0);
    if ~isequal(size(f0),[n,1])
        error('f must return a column vector with n rows.')
    end
    
    if isa(df,'function_handle')
        %df is already a function, great.
        g0 = df(x0);
        if ~isequal(size(g0),[n,n])
            error('g must return a matrix with size n by n.')
        end
    else
        error('Analytical gradient must be supplied.')
    end
    
    xVec = NaN(n,maxIt);
    fVec = NaN(n,maxIt);
    gVec = NaN(n,n,maxIt);
    deltaXvec = NaN(n,maxIt);
    
    xVec(:,1) = x0;
    fVec(:,1) = f0;
    gVec(:,:,1) = df(x0);
    
    i = 1;
    toContinue = true;
    while toContinue
        
        if norm(fVec(:,i)) < fTol 
            xf = xVec(:,i);
            ff = fVec(:,i);
            toContinue = false;
            toDirectionSearch = false;
            fprintf('done i = %i, f = %3.3e\n',i,norm(fVec(:,i)))
        else
            toDirectionSearch = true;
        end
        
        if toDirectionSearch
            gVec(:,:,i) = df(xVec(:,i));
            
            deltaXvec(:,i) = linsolve(gVec(:,:,i),-fVec(:,i));

            fCandidate = f(xVec(:,i) + deltaXvec(:,i));
            
            toLineSearch = norm(fCandidate) > norm(fVec(:,i));

            while toLineSearch
                deltaXvec(:,i) = deltaXvec(:,i)/2;

                fCandidate = f(xVec(:,i) + deltaXvec(:,i));

                if norm(deltaXvec(:,i)) < xTol
                    toContinue = false;
                    xf = xVec(:,i);
                    ff = fVec(:,i);
                    break
                end

                toLineSearch = norm(fCandidate) > norm(fVec(:,i));
            end

            fVec(:,i+1) = fCandidate;
            xVec(:,i+1) = xVec(:,i) + deltaXvec(:,i);

            i = i + 1;
            if i > maxIt
                toContinue = false;
                xf = NaN;
                ff = NaN;
                keyboard
            end
        end
        
    end
    
    if strcmp(class(figHandle),'matlab.ui.Figure') && n == 2
        customColorMap = parula(i);
        
        figure(figHandle)
        hold on
        for j=1:min(i,maxIt)
            pVec = [xVec(:,j),xVec(:,j) + deltaXvec(:,j)]
            plot(pVec(1,:),pVec(2,:),'k','linewidth',0.25)
        end
        scatter3(xVec(1,1:i),xVec(2,1:i),sqrt(sum(fVec(:,1:i).^2)),50,customColorMap,'filled')
        
        if all(~isnan(xf))
            title(sprintf('Solution at (x,y) = (%3.3f,%3.3f) found after %i loops.',xf(1),xf(2),i))
        else
            title('Solution did not converge')
        end
        
    end
    
end


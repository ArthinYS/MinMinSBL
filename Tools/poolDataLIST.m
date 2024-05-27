function yout = poolDataLIST(yin,ahat,nVars,polyorder,usesine)

% Copyright 2015, All Rights Reserved
% Code by Steven L. Brunton
% For Paper, "Discovering Governing Equations from Data: 
%        Sparse Identification of Nonlinear Dynamical Systems"
% by S. L. Brunton, J. L. Proctor, and J. N. Kutz

n = size(yin,1);

ind = 1;
% poly order 0
yout{ind,1} = ['1'];
ind = ind+1;

% poly order 1
for i=1:nVars
    yout{ind,1} =[ yin{i}];
    ind = ind+1;
end

if(polyorder>=2)
    % poly order 2
    for i=1:nVars
        for j=i:nVars
            yout{ind,1} = [yin{i},yin{j}];
            ind = ind+1;
        end
    end
end

if(polyorder>=3)
    % poly order 3
    for i=1:nVars
        for j=i:nVars
            for k=j:nVars
                yout{ind,1} = [yin{i},yin{j},yin{k}];
                ind = ind+1;
            end
        end
    end
end

if(polyorder>=4)
    % poly order 4
    for i=1:nVars
        for j=i:nVars
            for k=j:nVars
                for l=k:nVars
                    yout{ind,1} = [yin{i},yin{j},yin{k},yin{l}];
                    ind = ind+1;
                end
            end
        end
    end
end

if(polyorder>=5)
    % poly order 5
    for i=1:nVars
        for j=i:nVars
            for k=j:nVars
                for l=k:nVars
                    for m=l:nVars
                        yout{ind,1} = [yin{i},yin{j},yin{k},yin{l},yin{m}];
                        ind = ind+1;
                    end
                end
            end
        end
    end
end

if(usesine)
    for k=1
        yout{ind,1} = ['sin(','x)'];
        ind = ind + 1;
          yout{ind,1} = ['sin(','y)'];
        ind = ind + 1;
          yout{ind,1} = ['sin(','z)'];
        ind = ind + 1;
        yout{ind,1} = ['cos(','x)'];
        ind = ind + 1;
        yout{ind,1} = ['cos(','y)'];
        ind = ind + 1;
        yout{ind,1} = ['cos(','z)'];
        ind = ind + 1;
    end
    for k=2:10;
        yout{ind,1} = ['sin(',num2str(k),'x)'];
        ind = ind + 1;
          yout{ind,1} = ['sin(',num2str(k),'y)'];
        ind = ind + 1;
          yout{ind,1} = ['sin(',num2str(k),'z)'];
        ind = ind + 1;
        yout{ind,1} = ['cos(',num2str(k),'x)'];
        ind = ind + 1;
        yout{ind,1} = ['cos(',num2str(k),'y)'];
        ind = ind + 1;
        yout{ind,1} = ['cos(',num2str(k),'z)'];
        ind = ind + 1;
    end
end


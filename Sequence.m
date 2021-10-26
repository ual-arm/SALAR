function h2 = Sequence(th2i)
    N=length(th2i);
    idx=zeros(N,1);
    resta=zeros(N-2,1);
    [~,j] = min(th2i);
    for i=1:N
        if mod(j+i+N,N)==0
            idx(i)=N;
        else
            idx(i)=mod(j+i+N,N);
        end
    end
    for i=1:N-2
        resta(i)=th2i(idx(i+1))-th2i(idx(i));
    end
    if min(resta)<=0
        h2=1;
    else
        h2=0;
    end
end

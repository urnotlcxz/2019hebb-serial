clear input hiddenneuron output sumwx hiddenneuronnew startpoint output M I history number flag expression s1 s2 s3
OUTPUT=2; 
HIDDEN=9;
STEP=1;

w=3*randn(HIDDEN*HIDDEN,1);
hiddenneuron=ones(1,HIDDEN);

len=0;
for path=1:1000
    for s=1:STEP
        for i=1:HIDDEN
            sumwx=0;
            for j=1:HIDDEN
                sumwx=sumwx+hiddenneuron(j)*w((i-1)*HIDDEN+j);
            end
            hiddenneuronnew(i)=exp(-1*sumwx*sumwx);%cos(sumwx);%(sin(sumwx)+1)/2;%
        end
        output=hiddenneuron(1:OUTPUT);
        if(output(OUTPUT)==1)
            output(OUTPUT)=output(OUTPUT)-1.0E-100;
        end
        hiddenneuron=hiddenneuronnew;

    end
    history(path,:)=output;
    plot(history);
    pause(0.001);
    len=len+1;
end

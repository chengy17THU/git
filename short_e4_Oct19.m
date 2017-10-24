% launched on Oct 19th
%% settings
% fix Eon0 according to XPS characterization
% however Eoff0 is no more fixed at 3*Eon0
nbond=1000;
lratio=0.3;
absci=0.00001:0.00001:0.0027;% mm/s
T=[303.15 402.15 471.15];
load ww.mat
EonO=285;EonOH=532;EonF=686;% Eon for C=O, C-OH and CnF
t1i=0;t2i=0;t3i=0;
abc_array=zeros(48,13);
% in 1 row of abc_array: a-b-cs under t1 t2 t3 respectively(3+3+3),
% corresponding [kbond,ratio,l0,k](4)
%% generate the abc array
% ratio of number between O, OH and F(namely C=O, C-OH and CnF) is unknown
orderabc=1;
for a=7:1:10
    for b=1:1:4
        for c=1:1:4
            if [a/b,b/c]~=[1,1]
                abc(orderabc,:)=[a,b,c];% no preallocation because of using size(orderabc)
                orderabc=orderabc+1;
            end
        end
    end
end
%% loop, no nesting
j=1;flagi=0;
input_t1=zeros((8*8*7*15),4);input_t2=zeros((8*8*7*15),4);input_t3=zeros((8*8*7*15),4);
flag=zeros(1000,48,13);% store the final results
for kbond=500:100:1200 % 8 times
    for ratio=1.5:0.5:5.0 % 8 times
        % "Eoff/Eon"s are different between 3 types of bonds theoretically
        EoffO=ratio*EonO;EoffOH=ratio*EonOH;EoffF=ratio*EonF;
        for k=7000:1000:13000 % 7 times
            for l0=0.1:0.2:2.9 % 15 times
                t1i=0;t2i=0;t3i=0;
                for orderabc=1:1:size(abc) % 48 times
                    a=abc(orderabc,1);b=abc(orderabc,2);c=abc(orderabc,3);
                    % this time, no transverse comparison between difftemps
                    % 30 Celsius/303.15K
                    [pu41O,favg41O]=pucal2(a*nbond,kbond,lratio,EoffO,k,EonO,k,T(1),absci,l0,1);
                    [pu41OH,favg41OH]=pucal2(b*nbond,kbond,lratio,3*EoffOH,k,EonOH,k,T(1),absci,l0,1);
                    [pu41F,favg41F]=pucal2(c*nbond,kbond,lratio,3*EoffF,k,EonF,k,T(1),absci,l0,1);
                    favg41=favg41O/1e6+favg41OH/1e6+favg41F/1e6;
                    pu41=(pu41O>=0)&(pu41OH>=0)&(pu41F>=0);
                    if pu41
                        if (favg41(3)>favg41(21))&&(favg41(21)>favg41(50))&&(favg41(87)>favg41(256))
                            if (favg41(256)/favg41(3)<0.75)&&(favg41(256)/favg41(3)>0.45)
                                t1i=t1i+1;
                                abc_array(orderabc,1:3)=[a,b,c];
                            end
                        end
                    end
                    % 129 celsius/402.15K
                    [pu42O,favg42O]=pucal2(a*nbond,kbond,lratio,EoffO,k,EonO,k,T(2),absci,l0,1);
                    [pu42OH,favg42OH]=pucal2(b*nbond,kbond,lratio,EoffOH,k,EonOH,k,T(2),absci,l0,1);
                    [pu42F,favg42F]=pucal2(c*nbond,kbond,lratio,EoffF,k,EonF,k,T(2),absci,l0,1);
                    favg42=favg42O/1e6+favg42OH/1e6+favg42F/1e6;
                    pu42=(pu42O>=0)&(pu42OH>=0)&(pu42F>=0);
                    if pu42
                        if (favg42(3)<favg42(21))&&(favg42(21)<favg42(50))&&(favg42(87)<favg42(256))
                            if (favg42(141)/favg42(10)<2.5)&&(favg42(141)/favg42(10)>1.5)
                                % if((favg42(256)-favg42(50))<=1*(favg42(50)-favg42(3)))&&(favg42(256)-favg42(50))>=0.7*(favg42(50)-favg42(3))
                                t2i=t2i+1;
                                abc_array(orderabc,4:6)=[a,b,c];
                                % end
                            end
                        end
                    end
                    %198 celsius/471.15K
                    [pu43O,favg43O]=pucal2(a*nbond,kbond,lratio,EoffO,k,EonO,k,T(3),absci,l0,1);
                    [pu43OH,favg43OH]=pucal2(b*nbond,kbond,lratio,EoffOH,k,EonOH,k,T(3),absci,l0,1);
                    [pu43F,favg43F]=pucal2(c*nbond,kbond,lratio,EoffF,k,EonF,k,T(3),absci,l0,1);
                    favg43=favg43O/1e6+favg43OH/1e6+favg43F/1e6;
                    pu43=(pu43O>=0)&(pu43OH>=0)&(pu43F>=0);
                    if pu43
                        if (favg43(3)<favg43(21))&&(favg43(21)<favg43(50))&&(favg43(87)<favg43(256))
                            if ((favg43(256)-favg43(50))<=1.2*(favg43(50)-favg43(3)))&&(favg43(256)-favg43(50))>=0.7*(favg43(50)-favg43(3))
                                % guarantee its linearity
                                % if favg43(256)/favg43(28)>=1.6 %
                                % conflicts with the former condition
                                    t3i=t3i+1;
                                    abc_array(orderabc,7:9)=[a,b,c];
                                % end
                            end
                        end
                    end
                end
                if t1i
                    input_t1(j,:)=[kbond,ratio,l0,k];
                end
                if t2i
                    input_t2(j,:)=[kbond,ratio,l0,k];
                end
                if t3i
                    input_t3(j,:)=[kbond,ratio,l0,k];
                end
                if (t1i&&t2i)&&t3i 
                    % under this [kbond,ratio,l0,k] condition, if
                    % (t1i&&t2i)&&t3i~=0, store ALL corresponding a-b-cs
                    abc_array(1,10:13)=[kbond,ratio,l0,k];
                    flagi=flagi+1;
                    flag(flagi,:,:)=abc_array(:,:);
                end
                j
                j=j+1;
            end
        end
    end
end
save 'short_e4_Oct19.mat'
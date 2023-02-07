load kl2;
load ks;
% kl, gengood, genbad, genmix, impgood, impbad, impmix
% ks,1-gggc,2-iggc,3-gwgc,4-iwgc,5-gsgc,6-isgc,7-ggbc,8-igbc,9-gwbc,10-iwbc,11-gsbc,12-isbc,13-ggmc,14-igmc,15-gwmc,16-iwmc,17-gsmc,18-ismc
clear modks;
clear kld;
modks{1,1}=ks{1,1};
modks{2,1}=ks{3,1};
modks{3,1}=ks{5,1};
modks{4,1}=ks{7,1};
modks{5,1}=ks{9,1};
modks{6,1}=ks{11,1};
modks{7,1}=ks{13,1};
modks{8,1}=ks{15,1};
modks{9,1}=ks{17,1};
modks{10,1}=ks{2,1};
modks{11,1}=ks{4,1};
modks{12,1}=ks{6,1};
modks{13,1}=ks{8,1};
modks{14,1}=ks{10,1};
modks{15,1}=ks{12,1};
modks{16,1}=ks{14,1};
modks{17,1}=ks{16,1};
modks{18,1}=ks{18,1};

% KL comparison
% p is theory
% q is approximation
for i=1:9
    for l=1:6
        for j=1:3
            clear p;
            clear ss;
            clear q;
            clear d;

            q=kl2{l*2-1,i};
            p=modks{(l-1)*3+j,1};
            d=0;
            
            if(length(p)>length(q))
                q(length(q)+1:length(p))=q(length(q));
            end
            for kk=1:length(p)
                if(p(kk)==0)
                    ss(kk)=0;
                else
                    ss(kk)=p(kk)*log2(p(kk)/q(kk));
                end
                d=d+ss(kk);
            end
            kld{1,i}(l,j)=d;
            
            clear p;
            clear ss;
            clear q;
            clear d;

            p=kl2{l*2-1,i};
            q=modks{(l-1)*3+j,1};
            d=0;
            if(length(p)>length(q))
                q(length(q)+1:length(p))=q(length(q));
            end
            for kk=1:length(p)
                if(p(kk)==0)
                    ss(kk)=0;
                else
                    ss(kk)=p(kk)*log2(p(kk)/q(kk));
                end
                d=d+ss(kk);
            end
            kld{2,i}(l,j)=d;
        end
    end
end

for i=1:9
    for j=1:6
        for k=1:3
            
            if(kld{2,i}(j,k) < kld{1,i}(j,k))
                kld{3,i}(j,k)=kld{2,i}(j,k);
                kld{4,i}(j,k)=1;
            else
                kld{3,i}(j,k)=kld{1,i}(j,k);
                kld{4,i}(j,k)=0;
            end
            
            if(kld{3,i}(j,k)<0)
                kld{3,i}(j,k)=max(kld{2,i}(j,k),kld{1,i}(j,k));
                kld{4,i}(j,k)=1-kld{4,i}(j,k);
            end
        end
    end
end

%euclidean
% kl, gengood, genbad, genmix, impgood, impbad, impmix
% ks,1-gggc,2-iggc,3-gwgc,4-iwgc,5-gsgc,6-isgc,7-ggbc,8-igbc,9-gwbc,10-iwbc,11-gsbc,12-isbc,13-ggmc,14-igmc,15-gwmc,16-iwmc,17-gsmc,18-ismc
clear dist
for i=1:9
    for l=1:6
        for j=1:3
            
%-------------------------euclidean---------------------------
            clear ss;
            q=kl2{l*2-1,i};
            p=modks{(l-1)*3+j,1};
            d=0;
            
            if(length(p)>length(q))
                q(length(q)+1:length(p))=q(length(q));
            end
            for kk=1:length(p)
                ss(kk)=q(kk)-p(kk);
                d=d+(ss(kk))*(ss(kk));
            end
            dist{1,i}(l,j)=sqrt(d);

%-------------------------cityblock---------------------------
            clear ss;
            q=kl2{l*2-1,i};
            p=modks{(l-1)*3+j,1};
            d=0;
            
            if(length(p)>length(q))
                q(length(q)+1:length(p))=q(length(q));
            end
            for kk=1:length(p)
                ss(kk)=q(kk)-p(kk);
                d=d+abs(ss(kk));
            end
            dist{2,i}(l,j)=d;
            
%-------------------------minkowski---------------------------
            clear ss;
            q=kl2{l*2-1,i};
            p=modks{(l-1)*3+j,1};
            d=0;
            tp=3;
            if(length(p)>length(q))
                q(length(q)+1:length(p))=q(length(q));
            end
            for kk=1:length(p)
                ss(kk)=abs(q(kk)-p(kk));
                d=d+power(ss(kk),tp);
            end
            dist{3,i}(l,j)=nthroot(d,tp);
%-------------------------chebyshev---------------------------
            clear ss;
            q=kl2{l*2-1,i};
            p=modks{(l-1)*3+j,1};
            d=0;
            
            if(length(p)>length(q))
                q(length(q)+1:length(p))=q(length(q));
            end
            for kk=1:length(p)
                ss(kk)=abs(q(kk)-p(kk));
            end
            dist{4,i}(l,j)=max(ss);
            
            %-------------------------sorensen---------------------------
            clear ss;
            clear ss2;
            q=kl2{l*2-1,i};
            p=modks{(l-1)*3+j,1};
            d=0;
            d2=0;
            if(length(p)>length(q))
                q(length(q)+1:length(p))=q(length(q));
            end
            
            for kk=1:length(p)
                ss(kk)=abs(p(kk)-q(kk));
                ss2(kk)=abs(p(kk)+q(kk));
                
                d=d+ss(kk);
                d2=d2+ss2(kk);
            end
            dist{5,i}(l,j)=d/d2;
            clear ss2;
            clear d2;
 %-------------------------gower---------------------------
            clear ss;
            q=kl2{l*2-1,i};
            p=modks{(l-1)*3+j,1};
            d=0;
            
            if(length(p)>length(q))
                q(length(q)+1:length(p))=q(length(q));
            end
            
            for kk=1:length(p)
                ss(kk)=abs(p(kk)-q(kk));
                d=d+ss(kk);
            end
            dist{6,i}(l,j)=d/length(p);

 %-------------------------soergel---------------------------
            clear ss;
            clear ss2;
            q=kl2{l*2-1,i};
            p=modks{(l-1)*3+j,1};
            d=0;
            d2=0;
            
            if(length(p)>length(q))
                q(length(q)+1:length(p))=q(length(q));
            end
            
            for kk=1:length(p)
                ss(kk)=abs(p(kk)-q(kk));
                ss2(kk)=max(p(kk),q(kk));
                d=d+ss(kk);
                d2=d2+ss2(kk);
            end
            dist{7,i}(l,j)=d/d2;
            
 %-------------------------Kulczynski---------------------------
            clear ss;
            clear ss2;
            q=kl2{l*2-1,i};
            p=modks{(l-1)*3+j,1};
            d=0;
            d2=0;
            
            if(length(p)>length(q))
                q(length(q)+1:length(p))=q(length(q));
            end
            
            for kk=1:length(p)
                ss(kk)=abs(p(kk)-q(kk));
                ss2(kk)=min(p(kk),q(kk));
                d=d+ss(kk);
                d2=d2+ss2(kk);
            end
            dist{8,i}(l,j)=d/d2;
            
 %-------------------------canberra---------------------------
            clear ss;
            q=kl2{l*2-1,i};
            p=modks{(l-1)*3+j,1};
            d=0;
            
            if(length(p)>length(q))
                q(length(q)+1:length(p))=q(length(q));
            end
            
            for kk=1:length(p)
                ss(kk)=abs(p(kk)-q(kk))/(p(kk)+q(kk));
                if(abs(p(kk)-q(kk))== 0 && (p(kk)+q(kk))==0)
                    ss(kk)=0;
                end
                d=d+ss(kk);
            end
            dist{9,i}(l,j)=d;
        
 %-------------------------lorentzian---------------------------
            clear ss;
            q=kl2{l*2-1,i};
            p=modks{(l-1)*3+j,1};
            d=0;
            
            if(length(p)>length(q))
                q(length(q)+1:length(p))=q(length(q));
            end
            
            for kk=1:length(p)
                ss(kk)=log(1+abs(p(kk)-q(kk)));
                d=d+ss(kk);
            end
            dist{10,i}(l,j)=d;
            
                    
 %-------------------------intersection---------------------------
            clear ss;
            q=kl2{l*2-1,i};
            p=modks{(l-1)*3+j,1};
            d=0;
            
            if(length(p)>length(q))
                q(length(q)+1:length(p))=q(length(q));
            end
            
            for kk=1:length(p)
                ss(kk)=min(p(kk),q(kk));
                d=d+ss(kk);
            end
            dist{11,i}(l,j)=1-d;
            
             %-------------------------wave hedges---------------------------
            clear ss;
            q=kl2{l*2-1,i};
            p=modks{(l-1)*3+j,1};
            d=0;
            
            if(length(p)>length(q))
                q(length(q)+1:length(p))=q(length(q));
            end
            
            for kk=1:length(p)
                ss(kk)=1-(min(p(kk),q(kk))/max(p(kk),q(kk)));
                if(isnan(ss(kk)))
                    ss(kk)=0;
                end
                d=d+ss(kk);
            end
            dist{12,i}(l,j)=1-d;
            
             %-------------------------czekanowski---------------------------
            clear ss;
            clear ss2;
            q=kl2{l*2-1,i};
            p=modks{(l-1)*3+j,1};
            d=0;
            d2=0;
            if(length(p)>length(q))
                q(length(q)+1:length(p))=q(length(q));
            end
            
            for kk=1:length(p)
                ss(kk)=min(p(kk),q(kk));
                ss2(kk)=(p(kk)+q(kk));
                d=d+ss(kk);
                d2=d2+ss2(kk);
            end
            dist{13,i}(l,j)=1-(2*d/d2);
            
            clear ss2;
            clear d2;
            
                         %-------------------------motyka---------------------------
            clear ss;
            clear ss2;
            q=kl2{l*2-1,i};
            p=modks{(l-1)*3+j,1};
            d=0;
            d2=0;
            if(length(p)>length(q))
                q(length(q)+1:length(p))=q(length(q));
            end
            
            for kk=1:length(p)
                ss(kk)=min(p(kk),q(kk));
                ss2(kk)=(p(kk)+q(kk));
                d=d+ss(kk);
                d2=d2+ss2(kk);
            end
            dist{14,i}(l,j)=1-(d/d2);
            
            clear ss2;
            clear d2;
            
                         %-------------------------kulczynski---------------------------
            clear ss;
            clear ss2;
            q=kl2{l*2-1,i};
            p=modks{(l-1)*3+j,1};
            d=0;
            d2=0;
            if(length(p)>length(q))
                q(length(q)+1:length(p))=q(length(q));
            end
            
            for kk=1:length(p)
                ss(kk)=min(p(kk),q(kk));
                ss2(kk)=abs(p(kk)-q(kk));
                d=d+ss(kk);
                d2=d2+ss2(kk);
            end
            dist{15,i}(l,j)=d/d2;
            
            clear ss2;
            clear d2;
            
                         %-------------------------ruzicka---------------------------
            clear ss;
            clear ss2;
            q=kl2{l*2-1,i};
            p=modks{(l-1)*3+j,1};
            d=0;
            d2=0;
            if(length(p)>length(q))
                q(length(q)+1:length(p))=q(length(q));
            end
            
            for kk=1:length(p)
                ss(kk)=min(p(kk),q(kk));
                ss2(kk)=max(p(kk),q(kk));
                d=d+ss(kk);
                d2=d2+ss2(kk);
            end
            dist{16,i}(l,j)=d/d2;
            
            clear ss2;
            clear d2;
            
                         %-------------------------tanimoto---------------------------
            clear ss;
            clear ss2;
            clear ss3;
            q=kl2{l*2-1,i};
            p=modks{(l-1)*3+j,1};
            d=0;
            d2=0;
            d3=0;
            if(length(p)>length(q))
                q(length(q)+1:length(p))=q(length(q));
            end
            
            for kk=1:length(p)
%                 ss(kk)=max(p(kk),q(kk))-min(p(kk),q(kk));
%                 ss2(kk)=max(p(kk),q(kk));

                  ss(kk)=p(kk);
                  ss2(kk)=q(kk);
                  ss3(kk)=min(p(kk),q(kk));
                  
                  d=d+ss(kk);
                  d2=d2+ss2(kk);
                  d3=d3+ss3(kk);
%                 d=d+ss(kk);
%                 d2=d2+ss2(kk);
            end
%             dist{17,i}(l,j)=d/d2;
            dist{17,i}(l,j)=(ss(kk)+ss2(kk)-2*ss3(kk))/(ss(kk)+ss2(kk)-ss3(kk));
            clear ss3;
            clear d3;
            clear ss2;
            clear d2;
            
                        
                         %-------------------------innerproduct---------------------------
            clear ss;
            q=kl2{l*2-1,i};
            p=modks{(l-1)*3+j,1};
            d=0;
            
            if(length(p)>length(q))
                q(length(q)+1:length(p))=q(length(q));
            end
            
            for kk=1:length(p)
                ss(kk)=p(kk)*q(kk);
                d=d+ss(kk);
            end
            dist{18,i}(l,j)=d;
            
                                     %-------------------------harmonic  mean---------------------------
            clear ss;
            q=kl2{l*2-1,i};
            p=modks{(l-1)*3+j,1};
            d=0;
            
            if(length(p)>length(q))
                q(length(q)+1:length(p))=q(length(q));
            end
            
            for kk=1:length(p)
                ss(kk)=(p(kk)*q(kk))/(p(kk)+q(kk));
                if(isnan(ss(kk)))
                    ss(kk)=0;
                end
                d=d+ss(kk);
            end
            dist{19,i}(l,j)=2*d;
            
                                     %-------------------------cosine---------------------------
            clear ss;
            clear ss2;
            clear ss3;
            clear d2;
            clear d3;
            q=kl2{l*2-1,i};
            p=modks{(l-1)*3+j,1};
            d=0;
            d2=0;
            d3=0;
            if(length(p)>length(q))
                q(length(q)+1:length(p))=q(length(q));
            end
            
            for kk=1:length(p)
                ss(kk)=p(kk)*q(kk);
                ss2(kk)=p(kk)*p(kk);
                ss3(kk)=q(kk)*q(kk);
                d=d+ss(kk);
                d2=d2+ss2(kk);
                d3=d3+ss3(kk);
            end
            dist{20,i}(l,j)=d/(sqrt(d2)*sqrt(d3));
            clear d2;
            clear d3;
            clear ss2;
            clear ss3;
                        
             %-------------------------kumar-hassebrook---------------------------
            clear ss;
            clear ss2;
            clear ss3;
            clear d2;
            clear d3;
            q=kl2{l*2-1,i};
            p=modks{(l-1)*3+j,1};
            d=0;
            d2=0;
            d3=0;
            if(length(p)>length(q))
                q(length(q)+1:length(p))=q(length(q));
            end
            
            for kk=1:length(p)
                ss(kk)=p(kk)*q(kk);
                ss2(kk)=p(kk)*p(kk);
                ss3(kk)=q(kk)*q(kk);
                d=d+ss(kk);
                d2=d2+ss2(kk);
                d3=d3+ss3(kk);
            end
            dist{21,i}(l,j)=d/(d2+d3-d);
            clear d2;
            clear d3;
            clear ss2;
            clear ss3;
            
             %-------------------------jaccard---------------------------
            clear d;
            clear ss;
            clear ss2;
            clear ss3;
            clear d2;
            clear d3;
            q=kl2{l*2-1,i};
            p=modks{(l-1)*3+j,1};
            d=0;
            d2=0;
            d3=0;
            if(length(p)>length(q))
                q(length(q)+1:length(p))=q(length(q));
            end
            
            for kk=1:length(p)
                ss(kk)=p(kk)*q(kk);
                ss2(kk)=p(kk)*p(kk);
                ss3(kk)=q(kk)*q(kk);
                d=d+ss(kk);
                d2=d2+ss2(kk);
                d3=d3+ss3(kk);
            end
            dist{22,i}(l,j)=1-(d/(d2+d3-d));
            clear d2;
            clear d3;
            clear ss2;
            clear ss3;
            
             %-------------------------dice---------------------------
            clear ss;
            clear ss2;
            clear ss3;
            clear d2;
            clear d3;
            q=kl2{l*2-1,i};
            p=modks{(l-1)*3+j,1};
            d=0;
            d2=0;
            d3=0;
            if(length(p)>length(q))
                q(length(q)+1:length(p))=q(length(q));
            end
            
            for kk=1:length(p)
                ss(kk)=p(kk)*q(kk);
                ss2(kk)=p(kk)*p(kk);
                ss3(kk)=q(kk)*q(kk);
                d=d+ss(kk);
                d2=d2+ss2(kk);
                d3=d3+ss3(kk);
            end
            dist{23,i}(l,j)=1-((2*d)/(d2+d3));
            clear d2;
            clear d3;
            clear ss2;
            clear ss3;
                         %-------------------------fidelity---------------------------
            clear ss;
            q=kl2{l*2-1,i};
            p=modks{(l-1)*3+j,1};
            d=0;
            if(length(p)>length(q))
                q(length(q)+1:length(p))=q(length(q));
            end
            
            for kk=1:length(p)
                ss(kk)=sqrt(p(kk)*q(kk));
                d=d+ss(kk);
            end
            dist{24,i}(l,j)=d;
            
             %-------------------------Bhattacharyya---------------------------
            clear ss;
            q=kl2{l*2-1,i};
            p=modks{(l-1)*3+j,1};
            d=0;
            if(length(p)>length(q))
                q(length(q)+1:length(p))=q(length(q));
            end
            
            for kk=1:length(p)
                ss(kk)=sqrt(p(kk)*q(kk));
                d=d+ss(kk);
            end
            dist{25,i}(l,j)=-log(d);
            
             %-------------------------Hellinger---------------------------
            clear ss;
            q=kl2{l*2-1,i};
            p=modks{(l-1)*3+j,1};
            d=0;
            if(length(p)>length(q))
                q(length(q)+1:length(p))=q(length(q));
            end
            
            for kk=1:length(p)
                %ss(kk)=sqrt(p(kk)*q(kk));
                ss(kk)=(sqrt(p(kk))-sqrt(q(kk)))^2;
                d=d+ss(kk);
            end
           % dist{26,i}(l,j)=2*sqrt(1-d);
            dist{26,i}(l,j)=sqrt(2*d);
                         %-------------------------matusita---------------------------
            clear ss;
            q=kl2{l*2-1,i};
            p=modks{(l-1)*3+j,1};
            d=0;
            if(length(p)>length(q))
                q(length(q)+1:length(p))=q(length(q));
            end
            
            for kk=1:length(p)
                %ss(kk)=sqrt(p(kk)*q(kk));
                ss(kk)=(sqrt(p(kk))-sqrt(q(kk)))^2;
                d=d+ss(kk);
            end
            %dist{27,i}(l,j)=sqrt(2-2*d);
            dist{27,i}(l,j)=sqrt(d);
                        
                         %-------------------------squared-chord---------------------------
            clear ss;
            
            q=kl2{l*2-1,i};
            p=modks{(l-1)*3+j,1};
            d=0;
            if(length(p)>length(q))
                q(length(q)+1:length(p))=q(length(q));
            end
            
            for kk=1:length(p)
                ss(kk)=(sqrt(p(kk))-sqrt(q(kk)))^2;
                d=d+ss(kk);
            end
            dist{28,i}(l,j)=d;
            
                                     %-------------------------squared euclidean---------------------------
            clear ss;
            
            q=kl2{l*2-1,i};
            p=modks{(l-1)*3+j,1};
            d=0;
            if(length(p)>length(q))
                q(length(q)+1:length(p))=q(length(q));
            end
            
            for kk=1:length(p)
                ss(kk)=(p(kk)-q(kk))^2;
                d=d+ss(kk);
            end
            dist{29,i}(l,j)=d;
            
                                                 %-------------------------pearson---------------------------
            clear ss;
            
            q=kl2{l*2-1,i};
            p=modks{(l-1)*3+j,1};
            d=0;
            if(length(p)>length(q))
                q(length(q)+1:length(p))=q(length(q));
            end
            
            for kk=1:length(p)
                ss(kk)=((p(kk)-q(kk))^2)/q(kk);
                if(isnan(ss(kk)))
                    ss(kk)=0;
                end
                d=d+ss(kk);
            end
            dist{30,i}(l,j)=d;
            
                                                             %-------------------------neyman---------------------------
            clear ss;
            
            q=kl2{l*2-1,i};
            p=modks{(l-1)*3+j,1};
            d=0;
            if(length(p)>length(q))
                q(length(q)+1:length(p))=q(length(q));
            end
            
            for kk=1:length(p)
                ss(kk)=((p(kk)-q(kk))^2)/p(kk);
                if(isnan(ss(kk)))
                    ss(kk)=0;
                end
                d=d+ss(kk);
            end
            dist{31,i}(l,j)=d;
            
                        
         %-------------------------squared---------------------------
            clear ss;
            
            q=kl2{l*2-1,i};
            p=modks{(l-1)*3+j,1};
            d=0;
            if(length(p)>length(q))
                q(length(q)+1:length(p))=q(length(q));
            end
            
            for kk=1:length(p)
                ss(kk)=((p(kk)-q(kk))^2)/(p(kk)+q(kk));
                                if(isnan(ss(kk)))
                    ss(kk)=0;
                end
                d=d+ss(kk);
            end
            dist{32,i}(l,j)=d;
            
            %------------------------- Probabilistic Symmetric---------------------------
            clear ss;
            
            q=kl2{l*2-1,i};
            p=modks{(l-1)*3+j,1};
            d=0;
            if(length(p)>length(q))
                q(length(q)+1:length(p))=q(length(q));
            end
            
            for kk=1:length(p)
                ss(kk)=((p(kk)-q(kk))^2)/(p(kk)+q(kk));
                                if(isnan(ss(kk)))
                    ss(kk)=0;
                end
                d=d+ss(kk);
            end
            dist{33,i}(l,j)=2*d;
            
                        %-------------------------Divergence---------------------------
            clear ss;
            
            q=kl2{l*2-1,i};
            p=modks{(l-1)*3+j,1};
            d=0;
            if(length(p)>length(q))
                q(length(q)+1:length(p))=q(length(q));
            end
            
            for kk=1:length(p)
                ss(kk)=((p(kk)-q(kk))^2)/(p(kk)+q(kk))^2;
                                if(isnan(ss(kk)))
                    ss(kk)=0;
                end
                d=d+ss(kk);
            end
            dist{34,i}(l,j)=2*d;
            
                                    %-------------------------Clark---------------------------
            clear ss;
            
            q=kl2{l*2-1,i};
            p=modks{(l-1)*3+j,1};
            d=0;
            if(length(p)>length(q))
                q(length(q)+1:length(p))=q(length(q));
            end
            
            for kk=1:length(p)
                ss(kk)=((abs(p(kk)-q(kk)))/(p(kk)+q(kk)))^2;
                                if(isnan(ss(kk)))
                    ss(kk)=0;
                end
                d=d+ss(kk);
            end
            dist{35,i}(l,j)=sqrt(d);
            
                        
            %-------------------------additive symmetric---------------------------
            clear ss;
            
            q=kl2{l*2-1,i};
            p=modks{(l-1)*3+j,1};
            d=0;
            if(length(p)>length(q))
                q(length(q)+1:length(p))=q(length(q));
            end
            
            for kk=1:length(p)
                ss(kk)=((p(kk)-q(kk))^2*(p(kk)+q(kk)))/(p(kk)*q(kk));
                                if(isnan(ss(kk)))
                    ss(kk)=0;
                end
                d=d+ss(kk);
            end
            dist{36,i}(l,j)=d;
            
                        %-------------------------kullback-leibler---------------------------
            clear ss;
            
            q=kl2{l*2-1,i};
            p=modks{(l-1)*3+j,1};
            d=0;
            if(length(p)>length(q))
                q(length(q)+1:length(p))=q(length(q));
            end
            
            for kk=1:length(p)
                ss(kk)=p(kk)*log(p(kk)/q(kk));
                
                if(p(kk)==0)
                    ss(kk)=0;
                end
                d=d+ss(kk);
            end
            dist{37,i}(l,j)=d;
            
                                    %-------------------------Jeffreys---------------------------
            clear ss;
            
            q=kl2{l*2-1,i};
            p=modks{(l-1)*3+j,1};
            d=0;
            if(length(p)>length(q))
                q(length(q)+1:length(p))=q(length(q));
            end
            
            for kk=1:length(p)
                ss(kk)=(p(kk)-q(kk))*log(p(kk)/q(kk));
                if(isnan(ss(kk)))
                    ss(kk)=0;
                end
                d=d+ss(kk);
            end
            dist{38,i}(l,j)=d;
            
            
                                                %-------------------------K divergence---------------------------
            clear ss;
            
            q=kl2{l*2-1,i};
            p=modks{(l-1)*3+j,1};
            d=0;
            if(length(p)>length(q))
                q(length(q)+1:length(p))=q(length(q));
            end
            
            for kk=1:length(p)
                ss(kk)=p(kk)*log((2*p(kk))/(p(kk)+q(kk)));
                                if(isnan(ss(kk)))
                    ss(kk)=0;
                end
                d=d+ss(kk);
            end
            dist{39,i}(l,j)=d;
            
                                                            %-------------------------topsoe---------------------------
            clear ss;
            
            q=kl2{l*2-1,i};
            p=modks{(l-1)*3+j,1};
            d=0;
            if(length(p)>length(q))
                q(length(q)+1:length(p))=q(length(q));
            end
            
            for kk=1:length(p)
                ss(kk)=p(kk)*log((2*p(kk))/(p(kk)+q(kk))) + q(kk)*log((2*q(kk))/(p(kk)+q(kk))) ;
                                if(isnan(ss(kk)))
                    ss(kk)=0;
                end
                d=d+ss(kk);
            end
            dist{40,i}(l,j)=d;
            
                                                            %-------------------------jensen-shannon---------------------------
            clear ss;
            clear ss2;
            q=kl2{l*2-1,i};
            p=modks{(l-1)*3+j,1};
            d=0;
            d2=0;
            if(length(p)>length(q))
                q(length(q)+1:length(p))=q(length(q));
            end
            
            for kk=1:length(p)
                ss(kk)=p(kk)*log((2*p(kk))/(p(kk)+q(kk)));
                ss2(kk)=q(kk)*log((2*q(kk))/(p(kk)+q(kk)));
                
                if(isnan(ss2(kk)))
                    ss2(kk)=0;
                end

                if(isnan(ss(kk)))
                    ss(kk)=0;
                end
                d=d+ss(kk);
                d2=d2+ss2(kk);
            end
            dist{41,i}(l,j)=(1/2)*(d+d2);
            
                                                            %-------------------------jensen difference---------------------------
            clear ss;
            
            q=kl2{l*2-1,i};
            p=modks{(l-1)*3+j,1};
            d=0;
            if(length(p)>length(q))
                q(length(q)+1:length(p))=q(length(q));
            end
            
            for kk=1:length(p)
                ss(kk)=((p(kk)*log(p(kk))+q(kk)*log(q(kk)))/2)-(((p(kk)+q(kk))/2) * log((p(kk)+q(kk))/2));
                
                if(isnan(ss(kk)))
                    ss(kk)=0;
                end
                d=d+ss(kk);
            end
            dist{42,i}(l,j)=d;
            
                                                                        %-------------------------taneja---------------------------
            clear ss;
            
            q=kl2{l*2-1,i};
            p=modks{(l-1)*3+j,1};
            d=0;
            if(length(p)>length(q))
                q(length(q)+1:length(p))=q(length(q));
            end
            
            for kk=1:length(p)
                ss(kk)=((p(kk)+q(kk))/2)*(log((p(kk)+q(kk))/(2*sqrt(p(kk)*q(kk)))));
                
                if(isnan(ss(kk)))
                    ss(kk)=0;
                end
                d=d+ss(kk);
            end
            dist{43,i}(l,j)=d;
            
                                                                        %-------------------------kumar-johnson---------------------------
            clear ss;
            
            q=kl2{l*2-1,i};
            p=modks{(l-1)*3+j,1};
            d=0;
            if(length(p)>length(q))
                q(length(q)+1:length(p))=q(length(q));
            end
            
            for kk=1:length(p)
                ss(kk)=((p(kk)^2-q(kk)^2)^2)/(2*(p(kk)*q(kk))^(3/2));
                if(isnan(ss(kk)))
                    ss(kk)=0;
                end
                d=d+ss(kk);
            end
            dist{44,i}(l,j)=d;
            
                                                                        %-------------------------avg---------------------------
            clear ss;
            q=kl2{l*2-1,i};
            p=modks{(l-1)*3+j,1};
            d=0;
            if(length(p)>length(q))
                q(length(q)+1:length(p))=q(length(q));
            end
            
            for kk=1:length(p)
                ss(kk)=abs(p(kk)-q(kk));
                
                d=d+ss(kk);
            end
            dist{45,i}(l,j)=(d+max(ss))/2;
        end
    end
end
%Metric 1   subject1       goat    wolf    sheep        subject2
%           gen good
%           gen bad
%           gen mix
%           imp good
%           imp bad
%           imp mix
%Metric 2
%

clear output;
for i=1:size(dist,1)
    for j=1:size(dist,2)
        for k=1:size(dist{i,j},1)
            for l=1:size(dist{i,j},2)
                output((i-1)*size(dist{i,j},1)+k,(j-1)*size(dist{i,j},2)+l)=dist{i,j}(k,l);
            end
        end
    end
end

for i=1:45
    for j=1:3
        acc(i,j)=double(0);
    end
end
for i=1:size(output,1)
    for j=1:size(output,2)/3
        clear t;
        t(1)=output(i,(j-1)*3+1);
        t(2)=output(i,(j-1)*3+2);
        t(3)=output(i,(j-1)*3+3);
        
        [t1 t2]=min(t);
        outputTruth(i,(j-1)*3+t2)=1;
        
        if(j<=3 && t2==1)
            acc(ceil(i/6),1)=acc(ceil(i/6),1)+1;
        elseif(j>3&&j<=6 && t2==2)
            acc(ceil(i/6),2)=acc(ceil(i/6),2)+1;
        elseif(j>6&&j<=9 && t2==3)
            acc(ceil(i/6),3)=acc(ceil(i/6),3)+1;
        end
        
    end
end
for i=1:45
    acc(i,4)=sum(acc(i,1:3));
end
% 
% clear store;
% clear accd;
% w1=rand(45,1);
% w2=rand(45,1);
% stp=0.1;
% acc=0;
% cc=1;
% prev=0;
% while (acc<=0.60)
%     acc=0;
%     clear store;
%     for i=1:9  %should be 27
%         for q=1:6
%             yout1=0;
%             yout2=0;
%             ytar1=0;
%             ytar2=0;
% 
%             if(i<=3)
%                 ytar1=1;
%             elseif(i>3 && i<=6)
%                 ytar2=1;
%             end
% 
%             clear input;
%             for j=1:45
%                 input(j,1)=output((j-1)*6+q,i);
%             end
% 
%             for j=1:45
%                 yout1=yout1+w1(j)*output((j-1)*6+q,i);
%                 yout2=yout2+w2(j)*output((j-1)*6+q,i);
%             end
%             yy1=0;
%             yy2=0;
% 
%             if(yout1>0)
%                 yy1=1;
%             end
% 
%             if(yout2>0)
%                 yy2=1;
%             end
%             err1=ytar1-yy1;
%             err2=ytar2-yy2;
%             
%             wnew1=w1+stp*err1*input;
%             wnew2=w2+stp*err2*input;
% 
%             w1=wnew1;
%             w2=wnew2;
%         end
%     end
%     
%     for i=1:9
%         for q=1:6
%             yout1=0;
%             yout2=0;
%             ytar1=0;
%             ytar2=0;
% 
%             if(i<=3)
%                 ytar1=1;
%             elseif(i>3 && i<=6)
%                 ytar2=1;
%             end
% 
%             clear input;
%             for j=1:45
%                 input(j,1)=output((j-1)*6+q,i);
%             end
% 
%             for j=1:45
%                 yout1=yout1+w1(j)*output((j-1)*6+q,i);
%                 yout2=yout2+w2(j)*output((j-1)*6+q,i);
%             end
%             yy1=0;
%             yy2=0;
% 
%             if(yout1>0)
%                 yy1=1;
%             end
% 
%             if(yout2>0)
%                 yy2=1;
%             end
%             err1=ytar1-yy1;
%             err2=ytar2-yy2;
% 
%             if((err1==0) && (err2==0))
%                 acc=acc+1;
%                 
%                 store(i,q)=1;
%             end
% 
%         end
%     end
%     acc=acc/54;
%     if(acc> prev)
%        prev=acc
%     end
%    % accd(cc)=acc;
%    % cc=cc+1;
% end
outputt=[];
for i=1:45
    for j=1:6
        outputt((j-1)*45+i,1:27)=output((i-1)*6+j,1:27);
    end
end

TrainingSet=[];
GroupTrain=[];
k=0;
for j=1:6
    for i=1:9
        k=k+1;
        if (i<=3)
            GroupTrain=[GroupTrain; 1];
        elseif (i>3 & i<=6)
            
            GroupTrain=[GroupTrain; 2];
        else
            
            GroupTrain=[GroupTrain; 3];
        end
        
        TrainingSet(k,1:45)=outputt((j-1)*45+1:j*45,(i-1)*3+1);
        TrainingSet(k,46:45*2)=outputt((j-1)*45+1:j*45,(i-1)*3+2);
        TrainingSet(k,91:45*3)=outputt((j-1)*45+1:j*45,(i-1)*3+3);
    end
end
i=1;
j=1;
it=size(TrainingSet,1);
jt=size(TrainingSet,2);
while i<it
    while j<jt
        if(TrainingSet(i,j)==Inf)
            TrainingSet(:,j)=[];
            jt=jt-1;
            j=j-1;
        end
       
        j=j+1;
    end
    i=i+1;
end

for j=1:size(TrainingSet,1)
    results = multisvm(TrainingSet, GroupTrain, TrainingSet); 
    acc=0;
    for i=1:length(results)
        if(results(i)==GroupTrain(i))
            acc=acc+1;
        end
    end
    accd(j)=acc/length(results);
end
% TrainingSet=[ 1 10;2 20;3 30;4 40;5 50;6 66;3 30;4.1 42]; 
% TestSet=[3 34; 1 14; 2.2 25; 6.2 63]; 
% GroupTrain=[1;1;2;2;3;3;2;2]; 
% results = multisvm(TrainingSet, GroupTrain, TestSet); 
% disp('multi class problem'); 
% disp(results);




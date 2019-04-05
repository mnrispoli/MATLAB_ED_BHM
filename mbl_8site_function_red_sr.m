function [Ts,psiAllSave,Deig] = mbl_8site_function_red(Jx,U,W,DisordN,NT,isQP)


% experiment parameters
%phi_m=pi;
%phi_n=pi;

%{
Jx=0.038
J=Jx;
U=2.7.*J;
W=20.*J;
DisordN=199;
NT=151;
isQP=1;
%}


NSitesy=1;
NSitesx=8;
NSites=NSitesy*NSitesx;
NPart=9;
goldR=(1+sqrt(5))/2;
%goldR=2;
Aphi=0;
chi=-0;
qpdis=isQP;





xspec=round(NSitesx/2);
yspec=1;



HilbD=HilbDim(NSites,NPart);

basis = BasiskMake(HilbD,NSites,NPart);
%%



%{






Hjx=zeros(size(basis,1));
%Hjnx=Hjx;
%Hjnnx=Hjnx;

Hi=zeros(size(basis,1));

Hd=Hi;

jxterm=zeros([1 NSitesx]);
jxterm(1:2)=[1,-1];
jxterm=jxterm(toeplitz(1:numel(jxterm),[1 numel(jxterm):-1:2]));
%makes open boundary condition
jxterms=jxterm(2:end,:);
%jxterms=jxterm;
jxterms=vertcat(jxterms,fliplr(jxterms));


%other tunneling terms;
%{
jnxterm=zeros([1 NSitesx]);
jnxterm(1:3)=[1,0,-1];
jnxterm=jnxterm(toeplitz(1:numel(jnxterm),[1 numel(jnxterm):-1:2]));
%makes open boundary condition
jnxterms=jnxterm(3:end,:);
%jxterms=jxterm;
jnxterms=vertcat(jnxterms,fliplr(jnxterms));

jnnxterm=zeros([1 NSitesx]);
jnnxterm(1:4)=[1,0,0,-1];
jnnxterm=jnnxterm(toeplitz(1:numel(jnnxterm),[1 numel(jnnxterm):-1:2]));
%makes open boundary condition
jnnxterms=jnnxterm(4:end,:);
%jxterms=jxterm;
jnnxterms=vertcat(jnnxterms,fliplr(jnnxterms));



%%


%this works in particular for a ribbon geometry
%otherwise need to be more careful...
jyterm=zeros([1 NSites]);
jyterm(1)=1;
jyterm(1+NSitesx)=-1;
jyterm=jyterm(toeplitz(1:numel(jyterm),[1 numel(jyterm):-1:2]));
%jyterms=jyterm(2:end,:);
jyterms=jyterm;
%jyterms=vertcat(jyterms,fliplr(jyterms));
%}





tiltx=(NSitesx-1):-1:0;
tilt2dx=tiltx;
for kk=1:NSitesy-1
    tilt2dx=horzcat(tilt2dx,tiltx);
end
tiltx=fliplr(tiltx+1);

axDis=(1./goldR).*2.*pi.*tiltx;

tilty=(NSitesy-1):-1:0;
tilt2dy=tilty;
for kk=1:NSitesx-1
    tilt2dy=horzcat(tilt2dy,tilty);
end

%He(ii,jj)=sum(statei.*horzcat);

%{
d1s=[1 0 0 0 0 0 0 0];
d2s=[0 1 0 0 0 0 0 0];
d3s=[0 0 1 0 0 0 0 0];
d4s=[0 0 0 1 0 0 0 0];
d5s=[0 0 0 0 1 0 0 0];
d6s=[0 0 0 0 0 1 0 0];
d7s=[0 0 0 0 0 0 1 0];
d8s=[0 0 0 0 0 0 0 1];
%}

%axDis=[pi 0 0 0 0 0]
for ii=1:size(basis,1)
    ii
    statei=basis(ii,:);
    
    for jj=1:size(basis,1)
        statef=basis(jj,:);
        
        for kk=1:size(jxterms,1)
            %{
            yy=1;
            if kk>size(jxterms,1)/2
                yy=2;
            end
            %}
            stateix=statei;
            statefx=statef;
            %nearest
            if (stateix-statefx)==jxterms(kk,:)
                Hjx(ii,jj)=-sqrt(stateix(find(jxterms(kk,:)==1)))*sqrt(stateix(find(jxterms(kk,:)==-1))+1);
                [nii,njj]=find((statefx-stateix)==-1);
                [nii2,njj2]=find((statefx-stateix)==1);
                
                %{
                [bb,aa]=find(jxterms(kk,:)==1);
                nn=statefx(njj2);
 
                Hjx(ii,jj)=Hjx(ii,jj).*exp(-i*(chi*pi)*(nn)*sign(njj2-njj)).*exp(-i*(Aphi));
                Hjx(jj,ii)=conj(Hjx(ii,jj));
                %}
            end
        end
        
        %{
        %nextnearest
        for kk=1:size(jnxterms,1)
            %{
            yy=1;
            if kk>size(jxterms,1)/2
                yy=2;
            end
            %}
            stateix=statei;
            statefx=statef;
            %next nearest
            if (stateix-statefx)==jnxterms(kk,:)
                Hjnx(ii,jj)=-sqrt(stateix(find(jnxterms(kk,:)==1)))*sqrt(stateix(find(jnxterms(kk,:)==-1))+1);
                [nii,njj]=find((statefx-stateix)==-1);
                [nii2,njj2]=find((statefx-stateix)==1);
                
                [bb,aa]=find(jnxterms(kk,:)==1);
                nn=statefx(njj2);
 
                Hjnx(ii,jj)=Hjnx(ii,jj).*exp(-i*(chi*pi)*(nn)*sign(njj2-njj)).*exp(-i*(Aphi));
                Hjnx(jj,ii)=conj(Hjnx(ii,jj));
            end 
        end
        
        %nextnextnearest
        for kk=1:size(jnnxterms,1)
            %{
            yy=1;
            if kk>size(jxterms,1)/2
                yy=2;
            end
            %}
            stateix=statei;
            statefx=statef;
           %nextnextnearest 
           if (stateix-statefx)==jnnxterms(kk,:)
                Hjnnx(ii,jj)=-sqrt(stateix(find(jnnxterms(kk,:)==1)))*sqrt(stateix(find(jnnxterms(kk,:)==-1))+1);
                [nii,njj]=find((statefx-stateix)==-1);
                [nii2,njj2]=find((statefx-stateix)==1);
                
                [bb,aa]=find(jnnxterms(kk,:)==1);
                nn=statefx(njj2);
 
                Hjnnx(ii,jj)=Hjnnx(ii,jj).*exp(-i*(chi*pi)*(nn)*sign(njj2-njj)).*exp(-i*(Aphi));
                Hjnnx(jj,ii)=conj(Hjnnx(ii,jj));
            end 
            
            
        end
        
        %}

        
        if (statei-statef)==0
           Hi(ii,jj)=sum(statei.*(statei-1));

           
        end
        


                
    end
end



%depths=linspace(1,6,100)
%depths=

%}
%%
%for dd=1:length(depths)
%disordReals=rand(51,1)*1./1.6;
disordReals=rand(DisordN,1);



%% single particle mapping

mapss=cell(NPart,NSites);

for ss=1:NSites
    for nn=1:NPart+1
        mapss{nn,ss}=find(basis(:,ss)==nn-1);     
    end
end


%%

%load('8site_loadme_goldR.mat');
%load('8site_loadme.mat')

%load('8site_loadme_goldR2_pb.mat')
%load('8site_loadme.mat')
load('8s_9p_loadme.mat')

%load('Z:\DMD_templates\hor_fake_golden_lattice_17sites_X7_Y7_def-20_mode3_randpar0.15hologram\normalizedoffsets.mat')
%load('Z:\DMD_templates\hor_fake_golden_lattice_matt_temp_X7_Y7_def-20_mode3_randpar0.15hologram\normalizedoffsets.mat')

%load('data.mat')
%load('disord_pattern_DC2_def_n3.mat')
%load('Disord_UnNorm_Defocus_n3p3.mat');
%load('Disord_UnNorm_Defocus_n3p3.mat')
load('QPL_low_pxl_defocus_n30_DisordData.mat')


siteSave=savedData{3};
siteSave=siteSave./max(max(siteSave));

psiAllSave=cell(1,length(disordReals));

%%
for dd=1:length(disordReals)
    


%%

%init state
psiref=ones([1,NSitesx]);
psiref(1)=2;

psi_init=zeros([size(basis,1),1]);
psi11=find(ismember(basis,psiref,'rows')==1);
psi_init(psi11)=1;

Ts=logspace(-1.5,6,NT);


psi=zeros(length(Ts),size(basis,1));
%%
psi(1,:)=psi_init';
psiPerp=psi_init;
dt=Ts(2)-Ts(1);





theta=disordReals(dd).*2.*pi;



%%

AmpF=1;
AmpF=floor(200/DisordN);
%AmpF=1;

tryErr=0;
while tryErr==0
    try 

        if qpdis==1

            clear rvec
            dd
            Rvec=siteSave(8-(NSites/2-1):9+(NSites/2-1),floor(dd.*AmpF))';

            HdisRand=diag(sum((ones(size(basis,1),1)*Rvec).*basis,2));

            
            Ham=Jx.*Hjx+U./2.*Hi + ...  
                W.*HdisRand;
            
            theta=rand().*2.*pi;
        else
            %{
            dd
            HdisRand=diag(sum((ones(size(basis,1),1)*(2.*rand(1,NSitesx))).*basis,2));            
            Ham=Jx.*Hjx+U./2.*Hi+ ...  
                W.*HdisRand;
            %}
            
            clear rvec
            dd
             Rvec=siteSave(8-(NSites/2-1):9+(NSites/2-1),floor(dd.*AmpF))';
            
            szR=size(Rvec);
            
            Rvec=0.5.*cos([0:(NSites-1)]'.*goldR.*2.*pi+(2.*pi).*floor(dd.*AmpF)./DisordN)+0.5;
            Rvec=reshape(Rvec,szR);
            HdisRand=diag(sum((ones(size(basis,1),1)*Rvec).*basis,2));

            
            Ham=Jx.*Hjx+U./2.*Hi + ...  
                W.*HdisRand;
            
            theta=rand().*2.*pi;
            
   
        end


        [V,D]=eig(Ham);
        tryErr=1;
    catch
        tryErr=0;


    end
end
Deig(dd,:)=diag(D);
thetaStore(dd)=theta;

%%
%G.S. Init state



cns=conj(psi_init')*V;  

%%
for tt=1:length(Ts)

    psi(tt,:)=V*(cns'.*exp(diag(D).*i.*2.*pi.*Ts(tt)));


end
nbar=(abs(psi).^2)*basis;
nbar2=(abs(psi).^2)*(basis.^2);





psiAllSave{dd}=psi;

end




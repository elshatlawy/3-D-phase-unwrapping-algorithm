# 3-D phase unwrapping algorithm for Brain MRI images.

% (1) import the data - dataAy_MAG (magntiude images) + dataAy (phase images)
dataAy_MAG=dicom_read_tool
dataAy=dicom_read_tool

% (2) define the two phase echo-times te1 & te2
te1=2.8600
te2=6.7200

% (3) get a mask from the second echo-time magntiude images. (BET COMMAND)
BrainMask_te2 = load_nii('')

% get the size of all_data_img
all_data_img=double(dataAy.dataAy);
[lx ly lz]=size(all_data_img);
% get mask array of all_data_img. here
mask_data_img=BrainMask_te2.img;
idmask=find(BrainMask_te2.img(:));
mask_all_data_img=double(mask_data_img);
% reshapes the array 
% and change the range as phase is defined in the (-pi,pi) range
phase_data=(reshape(all_data_img,[lx,ly,lz/2,2])-2048)/2048*pi;
phase_data(:,:,:,1)=squeeze(phase_data(:,:,:,1)).*mask_all_data_img;
phase_data(:,:,:,2)=squeeze(phase_data(:,:,:,2)).*mask_all_data_img;
unWrappedPhase=zeros(size(phase_data));
 % squeeze the phase image (slice idS) first-echo-time to look at in any orientation 
phase_slc_eco1=squeeze(phase_data(:,:,:,1));
% squeeze the phase image (slice idS) second-echo-time to look at in any orientation
phase_slc_eco2=squeeze(phase_data(:,:,:,2));
% calculates the angle of the phase
phi211=angle(exp(j*phase_slc_eco2).*exp(-j*2.0*phase_slc_eco1));
% figure;imagesc(phi211);
% predict the phase images at the different TE using the unwrapping phase
predict_slc_eco1=phi211./(te2-2*te1)*te1;
wrapped_predict_slc_eco1=mod((predict_slc_eco1+pi),pi*2)-pi;
% unwrapping the phase map for different echo times.
wrappingThreshold=pi;
unwrapped_phase_slc_eco1=phase_slc_eco1;
for iw=1:fix(max(predict_slc_eco1(:)/2/pi))
    mask_pos_wrapped_slc_eco1=find(((predict_slc_eco1(:)-phase_slc_eco1(:))>((iw-1)*2*pi+wrappingThreshold))& ...
        ((predict_slc_eco1(:)-phase_slc_eco1(:))<=(iw*2*pi+wrappingThreshold)));
    unwrapped_phase_slc_eco1(mask_pos_wrapped_slc_eco1)=phase_slc_eco1(mask_pos_wrapped_slc_eco1)+iw*2*pi;
end
for iw=fix(min(predict_slc_eco1(:)/2/pi)):-1
    mask_nag_wrapped_slc_eco1=find(((predict_slc_eco1(:)-phase_slc_eco1(:))<((iw+1)*2*pi-wrappingThreshold))&...
        ((predict_slc_eco1(:)-phase_slc_eco1(:))>=(iw*2*pi-wrappingThreshold)));
    unwrapped_phase_slc_eco1(mask_nag_wrapped_slc_eco1)=phase_slc_eco1(mask_nag_wrapped_slc_eco1)+iw*2*pi;
end
unWrappedPhase(:,:,:,1)= unwrapped_phase_slc_eco1;
predict_slc_eco2=unwrapped_phase_slc_eco1/te1*te2;
wrapped_predict_slc_eco2=mod((predict_slc_eco2+pi),pi*2)-pi;
unwrapped_phase_slc_eco2=phase_slc_eco2;
for iw=1:fix(max(predict_slc_eco2(:)/2/pi))
    mask_pos_wrapped_slc_eco2=find(((predict_slc_eco2(:)-phase_slc_eco2(:))>((iw-1)*2*pi+wrappingThreshold))& ...
        ((predict_slc_eco2(:)-phase_slc_eco2(:))<=(iw*2*pi+wrappingThreshold)));
    unwrapped_phase_slc_eco2(mask_pos_wrapped_slc_eco2)=phase_slc_eco2(mask_pos_wrapped_slc_eco2)+iw*2*pi;
end
for iw=fix(min(predict_slc_eco2(:)/2/pi)):-1
    mask_nag_wrapped_slc_eco2=find(((predict_slc_eco2(:)-phase_slc_eco2(:))<((iw+1)*2*pi-wrappingThreshold))&...
        ((predict_slc_eco2(:)-phase_slc_eco2(:))>=(iw*2*pi-wrappingThreshold)));
    unwrapped_phase_slc_eco2(mask_nag_wrapped_slc_eco2)=phase_slc_eco2(mask_nag_wrapped_slc_eco2)+iw*2*pi;
end
% 
% figure;imagesc(unwrapped_phase_slc_eco2);
% figure;imagesc(phase_slc_eco2);
% figure;imagesc(mod((unwrapped_phase_slc_eco2+pi),pi*2)-pi);
unWrappedPhase(:,:,:,2)= unwrapped_phase_slc_eco2;


%%%%% Gradient Calculation %%%%

lzh=lz/2;
Gx_L=zeros([lx ly lzh]);
Gy_L=Gx_L;
Gz_L=Gx_L;
Gx_R=Gx_L;
Gy_R=Gx_L;
Gz_R=Gx_L;

 for i=2:lx     %when i=1  Gx_L = Gx_R
 Gx_L(i,:,:,1) = unWrappedPhase(i,:,:,1) - unWrappedPhase(i-1,:,:,1); 
 end
 
 for i=2:ly
 Gy_L(:,i,:,1) = unWrappedPhase(:,i,:,1) - unWrappedPhase(:,i-1,:,1);
 end
 
 for i=2:lzh
 Gz_L(:,:,i,1) = unWrappedPhase(:,:,i,1) - unWrappedPhase(:,:,i-1,1);
 end
 
 for i=1:lx-1 % when i=lx Gx_R = Gx_L
 Gx_R(i,:,:,1) = unWrappedPhase(i+1,:,:,1) - unWrappedPhase(i,:,:,1);
 end

 for i=1:ly-1
 Gy_R(:,i,:,1) = unWrappedPhase(:,i+1,:,1) - unWrappedPhase(:,i,:,1);
 end

 for i=1:lzh-1
 Gz_R(:,:,i,1) = unWrappedPhase(:,:,i+1,1) - unWrappedPhase(:,:,i,1);
 end
 
  Gx_L(1,:,:) = Gx_R(1,:,:);
  Gy_L(:,1,:) = Gy_R(:,1,:);
  Gz_L(:,:,1) = Gz_R(:,:,1);
  Gx_R(lx,:,:) = Gx_L(lx,:,:);
  Gy_R(:,ly,:) = Gy_L(:,ly,:);
  Gz_R(:,:,lzh) = Gz_L(:,:,lzh);
 

 Gabs=(sqrt(Gx_L.^2 + Gy_L.^2 + Gz_L.^2 + Gx_R.^2 + Gy_R.^2 + Gz_R.^2 ) ./ 2);
 
 
 %%% WEIGHTS %%%
 
  Magntiude = double(dataAy_MAG.dataAy(:,:,1:lzh));
  Magntiude_2 = double(dataAy_MAG.dataAy(:,:,lzh+1:lz));
  Mag_mean = sqrt(Magntiude(:,:,:) .* Magntiude_2(:,:,:));

Mx_L=zeros(size(Mag_mean));
My_L=zeros(size(Mag_mean));
Mz_L=zeros(size(Mag_mean));
Mx_R=zeros(size(Mag_mean));
My_R=zeros(size(Mag_mean));
Mz_R=zeros(size(Mag_mean));

    for m=2:lx
    Mx_L(m,:,:) = Magntiude_2(m-1,:,:);
    end
    
    for m=2:ly
    My_L(:,m,:) = Magntiude_2(:,m-1,:);
    end
    
    for m=2:lzh
    Mz_L(:,:,m) = Magntiude_2(:,:,m-1);
    end
    
    for m=1:lx-1
    Mx_R(m,:,:) = Magntiude_2(m+1,:,:);
    end
    
    for m=1:ly-1
    My_R(:,m,:) = Magntiude_2(:,m+1,:);
    end
    
    for m=1:lzh-1
    Mz_R(:,:,m) = Magntiude_2(:,:,m+1);
    end

    Mx_L(1,:,:) = Mx_R(1,:,:);
    My_L(:,1,:) = My_R(:,1,:);
    Mz_L(:,:,1) = Mz_R(:,:,1);
    Mx_R(lx,:,:) = Mx_L(lx,:,:);
    My_R(:,ly,:) = My_L(:,ly,:);
    Mz_R(:,:,lzh) = Mz_L(:,:,lzh);
  

%%%% FINAL GABS %%%%

 Gabs1 = (sqrt((Mx_L.*(Gx_L).^2))+((My_L.*(Gy_L).^2))+((Mz_L.*(Gz_L).^2))+((Mx_R.*(Gx_R).^2))+((My_R.*(Gy_R).^2))+((Mz_R.*(Gz_R).^2)));
 Gabs2 = (2 .* sqrt((Mx_L + My_L + Mz_L + Mx_R + My_R + Mz_R)));
 GabsF = (Gabs1 ./ Gabs2);
 
  
 sigma = std(Gabs(idmask(:)));
 Gradient_mask=Gabs<=sigma;
 NewFinal_Mask=(Gradient_mask(:,:,:) .* mask_all_data_img(:,:,:));
 mask_block=mask_all_data_img;
 for iz=1:lzh
     mask_block(:,:,iz)=mask_block(:,:,iz)*iz;
 end
 miniz=min(mask_block(:));
 maxiz=max(mask_block(:));
 mask_block=mask_all_data_img;
 for iy=1:ly
     mask_block(:,iy,:)=mask_block(:,iy,:)*iy;
 end
 miniy=min(mask_block(:));
 maxiy=max(mask_block(:));
 mask_block=mask_all_data_img;
 for ix=1:lx
     mask_block(ix,:,:)=mask_block(ix,:,:)*ix;
 end
 minix=min(mask_block(:));
 maxix=max(mask_block(:));
 mask_block=mask_all_data_img;
 xstep=4;
 mask_block(1:(minix-xstep),:,:)=0;
 mask_block((maxix+xstep):end,:,:)=0;
 mask_block(:,1:(miniy-xstep),:)=0;
 mask_block(:,(maxiy+xstep):end,:)=0; 
 mask_block(:,:,1:(miniz-xstep))=0;
 mask_block(:,:,(maxiz+xstep):end)=0;
 Gradient_mask = Gradient_mask.*mask_block;
   
 sigma_2 = std(GabsF(idmask(:)));
 Gradient_mask_2=GabsF<=sigma_2;
 NewFinal_Mask_2=(Gradient_mask_2(:,:,:) .* mask_all_data_img(:,:,:));
 Gradient_mask_2=Gradient_mask_2.*mask_block;
 
 %%%% Sorting analysis %%%%
 
[Gradient_mask_L,Gradient_mask_NUM]=bwlabeln(Gradient_mask,6);
LC_1=zeros(1,Gradient_mask_NUM+1);
for P=Gradient_mask_L(1:end); LC_1(P+1)=LC_1(P+1)+1; end
[LC_sorted, LC_index]=sort(-LC_1);
unWrappedPhase_Patched=unWrappedPhase(:,:,:,1);
unWrappedPhase_Patched(Gradient_mask_L~=LC_index(2)-1)=NaN;
unWrappedPhase_Patched(isnan(unWrappedPhase_Patched))=0;


[Gradient_mask_2_L,Gradient_mask_2_NUM]=bwlabeln(Gradient_mask_2,6);
unWrappedPhase_Patched_2=unWrappedPhase;
unWrappedPhase_Patched_2(Gradient_mask_2_L~=1)=NaN;
LC_1=zeros(1,Gradient_mask_2_NUM+1);
for P=Gradient_mask_2_L(1:end); LC_1(P+1)=LC_1(P+1)+1; end
[LC_sorted, LC_index]=sort(-LC_1);
unWrappedPhase_Patched_2=unWrappedPhase(:,:,:,1);
unWrappedPhase_Patched_2(Gradient_mask_2_L~=LC_index(2)-1)=NaN;
unWrappedPhase_Patched_2(isnan(unWrappedPhase_Patched_2))=0;
NewFinal_Mask_2(Gradient_mask_2_L~=LC_index(2)-1)=0;
 
 filter3d_mask=imgaussfilt3(NewFinal_Mask_2(:,:,:),1);
 filterMask = filter3d_mask>=0.5;
 filterMask=NewFinal_Mask_2(:,:,:);
 
 unwrappedphase_masked_NaN=unWrappedPhase_Patched_2;
 unwrappedphase_masked_NaN(filterMask==0)=NaN;
  
 unwrappedphase_masked_NaN_m3(filterMask==0)=NaN;
 unwrappedphase_masked_NaN_m3=medfilt3(squeeze(unwrappedphase_masked_NaN(:,:,:,1)),[3,3,3],'replicate');


%%%% Median filter %%%%

A=(unwrappedphase_masked_NaN(:,:,:,1));
filtered_unwrappedphase_masked_NaN=zeros(size(A));
IndNaN=isnan(A);
[x y z] = size(A);
Med = [];

for i=2:x-1 
	for j=2:y-1
		for k=2:z-1
		
			Med(1) = A(i-1,j,k);
			Med(2) = A(i+1,j,k);
			Med(3) = A(i,j-1,k);
			Med(4) = A(i,j+1,k);
			Med(5) = A(i,j,k-1);
			Med(6) = A(i,j,k+1);
			Med(7) = A(i,j,k);
			
			NotNaNs = ~isnan(Med);
			indNotNaN= find(NotNaNs==1);
			filtered_unwrappedphase_masked_NaN(i,j,k) = median(Med(indNotNaN));        
		end           
	end        
end
filtered_unwrappedphase_masked_NaN(find(IndNaN(:)==1))=NaN;


diff=A-filtered_unwrappedphase_masked_NaN;
sigma_fieldmap = nanstd(diff(idmask(:)));
mean_fieldmap = nanmean(diff(idmask(:)));

final_fieldmap_NaN=filtered_unwrappedphase_masked_NaN;

idfuwp=find(abs(diff(idmask(:)) - mean_fieldmap) <= sigma_fieldmap);
final_fieldmap_NaN(idmask(idfuwp))=A(idmask(idfuwp));
diff_aN=diff;diff_aN(idmask(idfuwp))=NaN;

%unWrappedPhase_NaN=unWrappedPhase(:,:,:,1);
%unWrappedPhase_NaN(filterMask==0)=NaN;

%unwrappedPhase_NaN=unwrappedPhase(:,:,:);
%unwrappedPhase_NaN(filterMask==0)=NaN;

final_fieldmap=final_fieldmap_NaN;
final_fieldmap(isnan(final_fieldmap))=0;
filtered_final_fieldmap=filtered_unwrappedphase_masked_NaN;
filtered_final_fieldmap(isnan(filtered_final_fieldmap))=0;
filterMask=double(filterMask);

%figure;imagesc(final_fieldmap(:,:,57))

    %%% END %%%     
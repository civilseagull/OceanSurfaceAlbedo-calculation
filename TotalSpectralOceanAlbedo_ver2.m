
 %Script made by Philipp M. Golovchenko (PhD)
 %For any questions - please contact - civilseagull@gmail.com
 
 
 %Spectral Albedo Calculation for Sea-Surface.
 %based on apparent-optical properties (proposed by Morel et al.),including
 %CHL(a) concentrations.
 %You will need for Input options:
 % Global Radiation Data
 % Chl(a)-concentrations 
 % Wavelength
 % Roughness -following Z.Jin scheme
 % 

 %Data read from csv-format file with Ktp and Local time (LT)

   % cases=[0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 ];%time steps
   cases=[8 9 10 11 12 13 14 15 16 17 18];
    % cases=[8 9 10];
    
    dataread=csvread('Global_Rad_3_Dates.csv');
    
for casenumb=1:length(cases); 
    
    caseworkingon= cases(casenumb);
    value_location=find(ismember(dataread(:,1),cases(casenumb)));
    
    LT=dataread(value_location,1); %Local Time
    GR= dataread(value_location,2); %Global Radiation
   
    
   %Basic Parameters Input%
    Lon=100; %Longitude
    Lat=5;   %Latitude
    d=29; %3 January- 3 day
    TZ=8; %Timezone
 
   %Physical Parameters Input%
   
       O=0;    %Surface Roughness
       W=440;  %Wavelength
       n0=1.341; %Refraction index of sea water
    
   %Pure Water coeffcients - by Raymond C. Smith and Karen S. Baker (1981)
   
       Aw=0.0145; %TYPE accordingly to table- for Wavelengths 420 =0.0153, 440 = 0,0145,  460 = 0.0156  560 = 0,0708, 660-0.4
       Bw=0.0049; % for 420 nm - 0.0061, 0.0049 - 440 nnm, 460 nm - 0.0041 0.0018 - 560 nm, 660nm - 0.0008
       
   %Chl-a Concentration (mg/m^3)
       Ch=1.1;
       
   %From Morel and Gentili 2001 -table value
       e=2.71828; %e -number
       
   %SOLAR Angles Option (IN DEGREES)%
   
     Dec=23.45.*sind((360./365).*(d-81)); %Declination angle
     B=(360/365).*(d).*(-81); %Special time factor
     EOT=9.87.*sind(2.*B)-7.53.*cosd(B)-1.5.*sind(B); %Equation of time
     LSTM=15.*TZ; %Local solar time meridian
     TC=4.*(Lon-LSTM)+(EOT); 
     LST=LT+(TC./60);  %Local Solar TIme
     HRA=15.*(LST-12); %Hour Angle
     SOLEL=asind(sind(Dec).*sind(Lat)+cosd(Dec).*cosd(Lat).*cosd(HRA)); %Solar Elevation Angle
     
     
   %Solar Zenith Angle Calculation%
   
     %Mu=sind(Lat).*sind(Dec)+cosd(Lat).*cosd(Dec).*cosd(HRA);%cos(SZA)
     Mu=sind(SOLEL);% COS(SZA)
     SZA=acosd(Mu); %Solar Zenith Angle
     
    %Solar Flux fractions (Direct and Diffused fractions for Spectral Albedo, by Jin et al. 2011)
    %Diffuse fraction of PAR%
    %Spitters model, 1986%

   %GLOBAL RADIATION INPUT options%
     
     Gsc=1367; %Solar Constant
     Gon=Gsc*(1+0.033*cosd((360*d)/365));%The exact irradiation incident on a surface
     Gext=Gon*(cosd(Lat)*cosd(Dec)*cosd(HRA)+sind(Lat)*sind(Dec));%extraterrestrial Radiation
     Kt=GR/Gext; %clearness index
     
   %DIFFUSE FRACTION OF GLOBAL RADIATION%
     
     %ERBS model (1982)
     %Kt<0.22, 
     Kd=1-(0.09*Kt);
     % 0.22<Kt<0.8  ,
    % Kd=0.9511 - 0.1604*Kt + 4.39*((Kt)^2) - 16.64*(Kt^3)+ 12.34*(Kt^4);
 
     
   %PAR ESTIMATION%
      PAR=0.49.*GR;
     %Spitters model, 1986%
      
     PARdif=(1+0.3*((1-Kd)^2))*Kd*PAR;
     PARdifF=PARdif/PAR; %Diffuse PAR fraction
     PARdirF=1 - PARdifF;%Direct PAR fraction
     
   %BIO-Optical Properties Input Calculation%
     
    %Backscattering of biological pigment -proposed in Morel and Maritorena (2001)%
   
     Bbp= 0.416*(Ch.^0.766)*(0.002+(1/100)*(0.50-0.25*(log(Ch))*((W./550)^(0.5*(log(Ch)-0.3)))));
     
     Nu=(0.5*Bw)/(0.5*Bw+Bbp); %the ratio of backscattering by water molecules to total backscattering
     B=0.6270-0.2227*Nu-0.0513*(Nu.^2)+(0.2465*Nu-0.3119)*Mu; %function of seawater and biological pigment backscattering
     
    %Absorption  coefficient of biological pigments (chl_a)
     
     ap1=0.18; %TYPE from 0.18 - 0.01 from Bricaud 1998
     ap= Ch.*ap1;
     a=(Aw + 0.06.*ap.*((Ch).^0.65)).*((1+0.2.*exp(-0.014.*(W-440))));%the absorption of Chl_a in m^-1
     Abp=0.06*a*(Ch.^0.65)+0.2*(0.00635+0.06*(Ch^0.65))*(e.^0.014*(440-W));
   
     %REFLECTANCE
     %represents an apparent optical property of seawater
     
       R0= B*((0.5*Bw+Bbp)/(Aw+Abp));
       
     %Fresnel reflectance (V.I. Haltrin, 2001)%
     Rf1= (cosd(SZA)-sqrt(n0.^2-sind(SZA).^2)/(cosd(SZA)+sqrt(n0.^2-sind(SZA).^2)).^2);
     Rf2= ((n0.^2).*cosd(SZA)-sqrt(n0.^2-sind(SZA).^2)/(n0.^2).*cosd(SZA)+sqrt(n0.^2-sind(SZA).^2)).^2;
     Rf=0.5.*(Rf1+Rf2);
     
     %Water-to-water reflectance at the air-water interface for upwelling diffuse incidence from water below
     rw=0.4817-0.0149.*O-0.207.*((O).^2.);
     
     
     %ALBEDO CALCULATION!!!%
     %Albedo divided following Jin (2011)%
    
                    %Variable names
                %ASdir -surface direct
                %ASdif -surface diffuse
                %AWdir - direct water
                %AWdif - duffuse water
     
     ASdir=Rf;
     %ASdif= -0.1482-0.012*O+0.1608*n0 - 0.0244*n0*O; %Isotropic incidence
     ASdif= -0.1479+0.1502*n0-0.0176*n0*O; %diffuse surface albedo under cloudy skies%
     AWdir= (R0*(1-rw)*(1-ASdir))/(1-rw*R0);
     %NOTE that AWdif=AWdir(but with Mu=0,676)!!! 
     
     B1=0.6270-0.2227*Nu-0.0513*(Nu.^2)+(0.2465*Nu-0.3119)*0.676;%with Mu=0,676
     R01= B1*((0.5*Bw+Bbp)/(Aw+Abp));  %with Mu=0,676
     
     AWdif=(R01*(1-rw)*(1-ASdir))/(1-rw*R01);
     
     
     %Direct ALBEDO component
      ALdir=ASdir+AWdir;
     %Diffuse ALBEDO component
      ALdif=ASdif+AWdif;
      
    %TOTAL SPECTRAL ALBEDO in PAR range for ocean surface%
     
      ALBEDO= (ALdir.*PARdirF +ALdif.*PARdifF);
      
    %Briegleb- 1986 scheme of Broadband OSA%

      AlBrieg= 0.026/(1.1*(Mu^(1.7))+0.065)+ 0.15*(Mu-0.1)*(Mu-0.5)*(Mu-1);

     
    %Output OPTIONS switcher
     
      outdata_TAlbedo(casenumb,:)= ALBEDO;
      outdata_KT(casenumb,:)= Kt;
      outdata_MU(casenumb,:)= Mu;
     
      outdataSZA(casenumb,:)=SZA;
      
      
    
     
end 
     
     
     
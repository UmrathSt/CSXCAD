function CSX = SetPBCExcitationWeight(CSX, name, sinweight, cosweight)
% function CSX = SetPBCExcitationWeight(CSX, name, sinweight, cosweight)
%
% Define PBC weighting functions for x-, y- and z-direction of excitation
%
% The functions can use the variables:
%   x,y,z
%   rho for the distance to z-axis
%   r   for the distance to origin
%   a   for alpha (as in cylindircal and spherical coord systems)
%   t   for theta (as in the spherical coord system
%   
%   all these variables are not weighted with the drawing unit defined by
%   the grid
% 
% example:
%     start=[0 0 0];
%     stop=[width height 0];
%     CSX = AddPBCExcitation(CSX,'excite',0,[1 1 0]);
%     sinweight{1} = '2*cos(0.0031416*x)*sin(0.0062832*y)';
%     sinweight{2} = '1*sin(0.0031416*x)*cos(0.0062832*y)';
%     sinweight{3} = 0;
%     cosweight{1} = '2*cos(0.0031416*x)*sin(0.0062832*y)';
%     cosweight{2} = '1*sin(0.0031416*x)*cos(0.0062832*y)';
%     cosweight{3} = 0;
%     CSX = SetPBCExcitationWeight(CSX,'excite',sinweight, cosweight);
%     CSX = AddBox(CSX,'excite',0 ,start,stop);
%
% See also AddExcitation, InitCSX, DefineRectGrid
% 
% CSXCAD matlab interface
% -----------------------
% author: Thorsten Liebig

if ~isfield(CSX,'Properties')
    error('CSXCAD::SetPBCExcitationWeight: no properties not found');
end
if ~isfield(CSX.Properties,'PBCExcitation')
    error('CSXCAD::SetPBCExcitationWeight: no excitation properties found');
end

pos=0;
for n=1:numel(CSX.Properties.PBCExcitation)
   if  strcmp(CSX.Properties.PBCExcitation{n}.ATTRIBUTE.Name, name)
       pos=n;
   end
end

if (pos==0)
    error('CSXCAD::SetPBCExcitationWeight: property not found');
    return;
end
CSX.Properties.PBCExcitation{pos}.SINWeight.ATTRIBUTE.X = sinweight{1};
CSX.Properties.PBCExcitation{pos}.SINWeight.ATTRIBUTE.Y = sinweight{2};
CSX.Properties.PBCExcitation{pos}.SINWeight.ATTRIBUTE.Z = sinweight{3};
CSX.Properties.PBCExcitation{pos}.COSWeight.ATTRIBUTE.X = cosweight{1};
CSX.Properties.PBCExcitation{pos}.COSWeight.ATTRIBUTE.Y = cosweight{2};
CSX.Properties.PBCExcitation{pos}.COSWeight.ATTRIBUTE.Z = cosweight{3};

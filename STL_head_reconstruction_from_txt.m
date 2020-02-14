%close all
clear all

figure(1)
clf
disp('Dock Figure!')
%pause(10)


scaling_factor = 0.5; %fill in PixelSpacing attribute from dicominfo of original file
model_base = 1; % adjust this to set the thickness of your model base (1mm is good)
model_z_offset = 0; %adjust this to set lowest point of model (where it merges with base)

iaz= [124 149 211]; %inter aureal zero in scan coordinates (only for chamber, ignore otherwise)

file = 'Y:\Projects\STS\Cornelius chamber planning\Cornelius MRIcron head model\rCO_20140423_STEREO_neurological_filled_and_smoothed5_indentsRL_center_corrected.txt';
o_file = 'Y:\Projects\STS\Cornelius chamber planning\Cornelius MRIcron head model\rCO_20140423_STEREO_neurological_filled_and_smoothed5_indentsRL_center_corrected.stl';

%following is optional, assign empty vector [] when not needed
voi_file = [];
voi_o_file= [];

% set to one if chamber construction (voi is required for position!)
make_chamber = 0;

if make_chamber == 1
chamber_height = 18;
chamber_width = 25;
chamber_width_inside = 19.1;
chamber_in_skull = 0.8;
chamber_offset_x = 4;
chamber_offset_y = 1;%2;
chamber_leg_radius = 15;
chamber_leg_thickness = 1.5;
%% 35deg chamber ? 
%% 0deg chamber ?
%% 55deg chamber ?
chamber_angle = -35 *pi/180;%-55*pi/180;%0;%-35 *pi/180;

ftk = 10; %frequencies to keep for bottom smoothing

c_file = 'sunny_chamber.stl';
c_file2 = 'sunny_chamber2.stl';
end



%% begin

disp('Reading file..')

X = textread(file);

X(:,1) = ((2+X(:,1))-min(X(:,1)))*scaling_factor;
X(:,2) = ((2+X(:,2))-min(X(:,2)))*scaling_factor;
X(:,3) = ((2+X(:,3))-min(X(:,3)))*scaling_factor;
X(:,3) = X(:,3) + model_z_offset;

diffs = abs(diff(X));
step = min(diffs(diffs~=0));
clear diffs

X2 = zeros(ceil(max(X(:,2)))+1,ceil(max(X(:,1)))+1);
for x=step:step:ceil(max(X(:,2)))+1
    p=0;
    for y=step:step:ceil(max(X(:,1)))+1
        p=p+1;
        if p==10
            fprintf('%s','.'); p=0;
        end
        index = find(abs(X(:,2)-x) < 1e-6 & abs(X(:,1)-y) <1e-6 );
        if isempty(index)
            X2(x/step,y/step) = model_base;
        else
            if max(X(index,3)) < model_base;
                X2(x/step,y/step) = model_base; 
            else
                X2(x/step,y/step) = max(X(index,3));
            end
        end
    end
    fprintf('%s\n',[num2str( (x/(ceil(max(X(:,2)))+1)*100),'%02.1f') '% '])
end
clear p

disp('Plot model')
hold off
surf(step:step:ceil(max(X(:,1)))+1,step:step:ceil(max(X(:,2)))+1,X2,'FaceAlpha',0.2)
axis equal


%return
disp('Writing STL ... please wait');

surf2stl(o_file,step,step,X2,'ascii');

clear x y index




% voi (optional)
if ~isempty(voi_file)
    
    disp('constructing volume of interest .. ');

O = textread(file);
Z = textread(voi_file);
o_file = voi_o_file;

Z(:,1) = ((2+Z(:,1))-min(O(:,1)))*scaling_factor;
Z(:,2) = ((2+Z(:,2))-min(O(:,2)))*scaling_factor;
Z(:,3) = ((2+Z(:,3))-min(O(:,3)))*scaling_factor;
Z(:,3) = Z(:,3) + model_z_offset;
clear O

Z = unique(Z,'rows');

T = delaunayn(Z,{'Qt','Qbb','Qc','Qz'});
% Limit circumradius of primitives
dt = TriRep(T,Z);
[~,rcc] = circumcenters(dt);
T = T(rcc < step*4,:);

F = [T(:,2:4); T(:,[1 3 4]); T(:,[1 2 4]); T(:,1:3)];
F = sort(F,2);
F = sortrows(F);
ii = find(~any(diff(F),2));
F([ii;ii+1],:) = [];

% Coordinates
x = Z(:,1);
y = Z(:,2);
z = Z(:,3);

% % Plot boundary faces
hold on
trisurf(F,x,y,z,'FaceColor','red','FaceAlpha',1);

vertices = [x y z];

v1=vertices(F(:,2),:)-vertices(F(:,1),:);
v2=vertices(F(:,3),:)-vertices(F(:,2),:);
% make the 3D normalized cross product
Norms = [(v1(:,2).*v2(:,3) - v1(:,3).*v2(:,2)) ...
         (v1(:,3).*v2(:,1) - v1(:,1).*v2(:,3)) ...
         (v1(:,1).*v2(:,2) - v1(:,2).*v2(:,1))];
       Norms_mag = sqrt(sum((Norms.*Norms)')');
       Norms(:,1) = Norms(:,1)./Norms_mag;
       Norms(:,2) = Norms(:,2)./Norms_mag;
       Norms(:,3) = Norms(:,3)./Norms_mag;
clear v1 v2 v3


v1(:,1:3)=vertices(F(:,1),1:3);
v2(:,1:3)=vertices(F(:,2),1:3);
v3(:,1:3)=vertices(F(:,3),1:3);

fid = fopen(o_file,'w');
fprintf(fid,'solid %s\n',o_file(1:length(o_file)-4));
%nf = length(F);
for k = 1:length(F)
fprintf(fid,'facet normal %5.5f %5.5f %5.5f\n outer loop\n vertex %5.5f %5.5f %5.5f\n vertex %5.5f %5.5f %5.5f\n vertex %5.5f %5.5f %5.5f\n endloop\n endfacet\n',...
     Norms(k,1),Norms(k,2),Norms(k,3), v1(k,1), v1(k,2), v1(k,3),v2(k,1), v2(k,2), v2(k,3),v3(k,1), v3(k,2), v3(k,3) );
end
fprintf(fid,'endsolid %s\n',o_file(1:length(o_file)-4));
fclose(fid);

clear Norms Norms_mag T dt k ii fid rcc v1 v2 v3 vertices

disp('done with voi ! ');

end


if make_chamber == 1

disp('pause')
pause(10)

disp('Chamber..')

center = [mean(x),mean(y),mean(z)];
center = center * (1/step);
center = round(center);
center = center / (1/step);

dist = zeros(size(X2));
for a=1:size(X2,1)
    for i=1:size(X2,2)
        xd=(i*step)-center(1);
        yd=(a*step)-center(2);
        zd= X2(a,i)-center(3);
        dist(a,i) = sqrt(xd*xd + yd*yd + zd*zd);
    end
end
        
for s=step:step:max(max(X2))
    itc(1) = center(1);
    itc(2) = center(2)+(sin(chamber_angle) * s);
    itc(3) = center(3)+(cos(chamber_angle) * s);
    getting_smaller = 0;
    for a=1:size(X2,1)
        for i=1:size(X2,2)
            xd=(i*step)-itc(1);
            yd=(a*step)-itc(2);
            zd= X2(a,i)-itc(3);
            distance = sqrt(xd*xd + yd*yd + zd*zd);
            if dist(a,i) > distance
                dist(a,i) = distance;
                getting_smaller = 1;
            end
        end
    end
    if getting_smaller == 0
        disp(['Found! ' int2str(s)])
        break;
    end
end
[~,index] = min(reshape(dist,numel(dist),1));
[a,i] = ind2sub(size(dist),index);

center_skull = [i*step a*step X2(a,i)] 

plot3([center(1),center_skull(1)],[center(2),center_skull(2)],[center(3),center_skull(3)],'r','LineWidth',5)
plot3(  [center_skull(1),center_skull(1)], ...
        [center_skull(2),center_skull(2)+(sin(chamber_angle) * chamber_height)], ...
        [center_skull(3),center_skull(3)+(cos(chamber_angle) * chamber_height)],'g','LineWidth',5)

oben = zeros(360*4 + 1,3);
oben(1,:) = [center_skull(1)+chamber_offset_x, ...
             center_skull(2)+(sin(chamber_angle) * chamber_height) + chamber_offset_y, ...
             center_skull(3)+(cos(chamber_angle) * chamber_height)];
%plot3(oben(1,1),oben(1,2),oben(1,3),'b*');
for i=1:360*4
    oben(i+1,:) = [ oben(1,1) + cos(i/4*pi/180)*chamber_width/2 , ...
                    oben(1,2) + sin(i/4*pi/180)*(chamber_width/2*cos(chamber_angle)), ...
                    oben(1,3) - sin(i/4*pi/180)*(chamber_width/2*sin(chamber_angle)) ];
    %plot3(oben(i+1,1),oben(i+1,2),oben(i+1,3),'k.');
end

drawnow

unten = zeros(360*4 + 1,3);
d=1;
for d=d:360*4 + 1
    for s=-chamber_height/step*1.5:0.1:0
        itc(1) = oben(d,1);
        itc(2) = oben(d,2)+(sin(chamber_angle) * s);
        itc(3) = oben(d,3)+(cos(chamber_angle) * s);
        getting_smaller = 1;
        for a=1:size(X2,1)
            for i=1:size(X2,2)
                xd=(i*step)-itc(1);
                yd=(a*step)-itc(2);
                zd= X2(a,i)-itc(3);
                distance = sqrt(xd*xd + yd*yd + zd*zd);
                if distance < chamber_in_skull
                    getting_smaller = 0;
                    break;
                end
            end
            if getting_smaller == 0
                break;
            end
        end
        if getting_smaller == 0
            disp(['Found by iteration: x=' num2str(itc(1)) ' y=' num2str(itc(2)) ' z=' num2str(itc(3)) ' (dist=' num2str(distance) ') for ' num2str((d-1)/4) '°'])
            unten(d,:) = [itc(1) itc(2) itc(3)] ;
            break;
        end
    end 
end
unten(1,:) = [unten(1,1), ...
             unten(1,2)+(sin(chamber_angle) * chamber_in_skull), ...
             unten(1,3)+(cos(chamber_angle) * chamber_in_skull)];
plot3(unten(:,1),unten(:,2),unten(:,3),'b.');
drawnow

oben_in = zeros(360*4 + 1,3);
oben_in(1,:) = oben(1,:);
for i=1:360*4
    oben_in(i+1,:) = [ oben_in(1,1) + cos(i/4*pi/180)*chamber_width_inside/2 , ...
                       oben_in(1,2) + sin(i/4*pi/180)*(chamber_width_inside/2*cos(chamber_angle)), ...
                       oben_in(1,3) - sin(i/4*pi/180)*(chamber_width_inside/2*sin(chamber_angle)) ];
end

unten_in = zeros(360*4 + 1,3);
d=1;
for d=d:360*4 + 1
    for s=-chamber_height/step*1.5:0.1:0
        itc(1) = oben_in(d,1);
        itc(2) = oben_in(d,2)+(sin(chamber_angle) * s);
        itc(3) = oben_in(d,3)+(cos(chamber_angle) * s);
        getting_smaller = 1;
        for a=1:size(X2,1)
            for i=1:size(X2,2)
                xd=(i*step)-itc(1);
                yd=(a*step)-itc(2);
                zd= X2(a,i)-itc(3);
                distance = sqrt(xd*xd + yd*yd + zd*zd);
                if distance < chamber_in_skull
                    getting_smaller = 0;
                    break;
                end
            end
            if getting_smaller == 0
                break;
            end
        end
        if getting_smaller == 0
            disp(['Found by iteration 2: x=' num2str(itc(1)) ' y=' num2str(itc(2)) ' z=' num2str(itc(3)) ' (dist=' num2str(distance) ') for ' num2str((d-1)/4) '°'])
            unten_in(d,:) = [itc(1) itc(2) itc(3)] ;
            break;
        end
    end 
end
unten_in(1,:) = [unten_in(1,1), ...
             unten_in(1,2)+(sin(chamber_angle) * chamber_in_skull), ...
             unten_in(1,3)+(cos(chamber_angle) * chamber_in_skull)];
plot3(unten_in(:,1),unten_in(:,2),unten_in(:,3),'b.');
drawnow

%crop
idx = unten(:,3)~=0 & unten_in(:,3)~=0;
oben = oben(logical(idx),:);
oben_in = oben_in(logical(idx),:);
paper_x = ([0.0:0.25:360]/360)*pi*chamber_width;
paper_x = paper_x(logical(idx));
unten = unten(logical(idx),:);
unten_in = unten_in(logical(idx),:);
%make chamber correct height
for i=1:size(oben,1)
    xd=oben_in(i,1)-unten_in(i,1);
    yd=oben_in(i,2)-unten_in(i,2);
    zd=oben_in(i,3)-unten_in(i,3);
    distance_in(i) = sqrt(xd*xd + yd*yd + zd*zd);
    xd=oben(i,1)-unten(i,1);
    yd=oben(i,2)-unten(i,2);
    zd=oben(i,3)-unten(i,3);
    distance_out(i) = sqrt(xd*xd + yd*yd + zd*zd);
end
oben_offset = max(chamber_height-min(distance_in),chamber_height-min(distance_out));

oben = [oben(:,1), ...
             oben(:,2)+(sin(chamber_angle) * oben_offset), ...
             oben(:,3)+(cos(chamber_angle) * oben_offset)];
oben_in = [oben_in(:,1), ...
             oben_in(:,2)+(sin(chamber_angle) * oben_offset), ...
             oben_in(:,3)+(cos(chamber_angle) * oben_offset)];
         
%plot oben

plot3(oben(1,1),oben(1,2),oben(1,3),'b*');
plot3([oben(1,1),unten(1,1)],[oben(1,2),unten(1,2)],[oben(1,3),unten(1,3)],'g-','LineWidth',3)
for i=2:size(oben,1)
    plot3(oben(i,1),oben(i,2),oben(i,3),'b.');
    plot3(oben_in(i,1),oben_in(i,2),oben_in(i,3),'b.');
    plot3([oben(i,1),oben_in(i,1)],[oben(i,2),oben_in(i,2)],[oben(i,3),oben_in(i,3)],'k:','LineWidth',0.1)
    plot3([oben(i,1),unten(i,1)],[oben(i,2),unten(i,2)],[oben(i,3),unten(i,3)],'k:','LineWidth',0.1)
end
drawnow


figure
hold off
plot(paper_x(2:size(paper_x,2)),distance_out(2:size(distance_out,2)))
hold on
plot(paper_x(2:size(paper_x,2)),distance_in(2:size(distance_in,2)),'b--')

O = textread(file);
from_iaz_x = (((unten(1,1)/scaling_factor)-2+min(O(:,1)))- iaz(1)) * scaling_factor
from_iaz_y = (((unten(1,2)/scaling_factor)-2+min(O(:,2)))- iaz(2)) * scaling_factor
from_iaz_z = (((unten(1,3)/scaling_factor)-2+min(O(:,3)))- iaz(3)) * scaling_factor;
from_iaz_z = from_iaz_z - model_z_offset

text(30,15,['Angle ' num2str(chamber_angle*180/pi) '°'],'FontSize',21)
text(30,10,['AP0 ' num2str(from_iaz_y)],'FontSize',20)
text(30,7,['ML ' num2str(from_iaz_x)],'FontSize',20)
text(30,4,['H ' num2str(from_iaz_z)],'FontSize',20)

mt_from_iaz_x = (((center(1)/scaling_factor)-2+min(O(:,1)))- iaz(1)) * scaling_factor
mt_from_iaz_y = (((center(2)/scaling_factor)-2+min(O(:,2)))- iaz(2)) * scaling_factor
mt_from_iaz_z = (((center(3)/scaling_factor)-2+min(O(:,3)))- iaz(3)) * scaling_factor;
mt_from_iaz_z = mt_from_iaz_z - model_z_offset

hold on
%fft smooting
ftc = floor(size(distance_in,2)/2) - ftk;
f = fft(distance_in);
f(floor(size(distance_in,2)/2)+1-ftc:floor(size(distance_in,2)/2)+ftc) = zeros(ftc*2,1);
smoothed_distance_in = real(ifft(f));
f = fft(distance_out);
f(floor(size(distance_out,2)/2)+1-ftc:floor(size(distance_out,2)/2)+ftc) = zeros(ftc*2,1);
smoothed_distance_out = real(ifft(f));

plot(paper_x(2:size(paper_x,2)),smoothed_distance_in(2:size(smoothed_distance_in,2)),'r--')
plot(paper_x(2:size(paper_x,2)),smoothed_distance_out(2:size(smoothed_distance_out,2)),'r')

axis equal
axis([0 max(paper_x) 0 ceil(max(max([distance_in distance_out])))])
drawnow


%% make first cad model
if ~all(unten(:,3)~=0 & unten_in(:,3)~=0)
    disp('Not all angles detected.. this is dangerous! Consider checking your model!')
end
unten_legs = zeros(360*4 + 1,3);
unten_legs = unten_legs(logical(idx),:);
for i=2:size(unten_legs,1)
    temp = unten(i,:)-unten_in(i,:);
    temp = temp/norm(temp);
    unten_legs(i,:) = (temp * chamber_leg_radius) + unten(i,:) + ...
        [0, sin(chamber_angle) * chamber_in_skull, cos(chamber_angle) * chamber_in_skull];
end

%oben
model = [1 0 0; 0 cos(chamber_angle) -sin(chamber_angle); 0 sin(chamber_angle) cos(chamber_angle)] ...
        * oben(1,:)';
%model = [model [1 0 0; 0 cos(chamber_angle) -sin(chamber_angle); 0 sin(chamber_angle) cos(chamber_angle)] ...
%        * oben(2:size(oben,1),:)' - repmat(model(:,1),1,size(oben,1)-1) ];
%unten
model = [model [1 0 0; 0 cos(chamber_angle) -sin(chamber_angle); 0 sin(chamber_angle) cos(chamber_angle)] ...
        * unten_in(2:size(unten_in,1),:)' - repmat(model(:,1),1,size(unten_in,1)-1) ];
model = [model [1 0 0; 0 cos(chamber_angle) -sin(chamber_angle); 0 sin(chamber_angle) cos(chamber_angle)] ...
        * unten(2:size(unten,1),:)' - repmat(model(:,1),1,size(unten,1)-1) ];
model = [model [1 0 0; 0 cos(chamber_angle) -sin(chamber_angle); 0 sin(chamber_angle) cos(chamber_angle)] ...
        * unten(2:size(unten,1),:)' - repmat(model(:,1)-[0;0;chamber_in_skull],1,size(unten,1)-1) ];
model = [model [1 0 0; 0 cos(chamber_angle) -sin(chamber_angle); 0 sin(chamber_angle) cos(chamber_angle)] ...
        * unten_legs(2:size(unten_legs,1),:)' - repmat(model(:,1),1,size(unten_legs,1)-1) ];

model(:,1) = [0;0;-chamber_height];

%rotate 180°
model = [1 0 0; 0 cos(pi) -sin(pi); 0 sin(pi) cos(pi)] ...
        * model(:,:) ;

model = model';

%fft smooth unten

ftc = floor(size(distance_in,2)/2) - ftk;

bis = size(model,1);
von = 1+bis - (size(unten_legs,1)-1);
f1 = fft(model(von:bis,1));
f2 = fft(model(von:bis,2));
f3 = fft(model(von:bis,3));
f1(floor(size(model(von:bis,1),1)/2)+1-ftc:floor(size(model(von:bis,1),1)/2)+ftc) = zeros(ftc*2,1);
f2(floor(size(model(von:bis,2),1)/2)+1-ftc:floor(size(model(von:bis,2),1)/2)+ftc) = zeros(ftc*2,1);
f3(floor(size(model(von:bis,3),1)/2)+1-ftc:floor(size(model(von:bis,3),1)/2)+ftc) = zeros(ftc*2,1);
model(von:bis,:) = [real(ifft(f1)) real(ifft(f2)) real(ifft(f3))];

bis1 = von-1;
von1 = 1+bis1 - (size(unten,1)-1);
f1 = fft(model(von1:bis1,1));
f2 = fft(model(von1:bis1,2));
f3 = fft(model(von1:bis1,3));
f1(floor(size(model(von1:bis1,1),1)/2)+1-ftc:floor(size(model(von1:bis1,1),1)/2)+ftc) = zeros(ftc*2,1);
f2(floor(size(model(von1:bis1,2),1)/2)+1-ftc:floor(size(model(von1:bis1,2),1)/2)+ftc) = zeros(ftc*2,1);
f3(floor(size(model(von1:bis1,3),1)/2)+1-ftc:floor(size(model(von1:bis1,3),1)/2)+ftc) = zeros(ftc*2,1);
model(von1:bis1,:) = [real(ifft(f1)) real(ifft(f2)) real(ifft(f3))];

bis2 = von1-1;
von2 = 1+bis2 - (size(unten,1)-1);
f1 = fft(model(von2:bis2,1));
f2 = fft(model(von2:bis2,2));
f3 = fft(model(von2:bis2,3));
f1(floor(size(model(von2:bis2,1),1)/2)+1-ftc:floor(size(model(von2:bis2,1),1)/2)+ftc) = zeros(ftc*2,1);
f2(floor(size(model(von2:bis2,2),1)/2)+1-ftc:floor(size(model(von2:bis2,2),1)/2)+ftc) = zeros(ftc*2,1);
f3(floor(size(model(von2:bis2,3),1)/2)+1-ftc:floor(size(model(von2:bis2,3),1)/2)+ftc) = zeros(ftc*2,1);
model(von2:bis2,:) = [real(ifft(f1)) real(ifft(f2)) real(ifft(f3))];

bis3 = von2-1;
von3 = 1+bis3 - (size(unten_in,1)-1);
f1 = fft(model(von3:bis3,1));
f2 = fft(model(von3:bis3,2));
f3 = fft(model(von3:bis3,3));
f1(floor(size(model(von3:bis3,1),1)/2)+1-ftc:floor(size(model(von3:bis3,1),1)/2)+ftc) = zeros(ftc*2,1);
f2(floor(size(model(von3:bis3,2),1)/2)+1-ftc:floor(size(model(von3:bis3,2),1)/2)+ftc) = zeros(ftc*2,1);
f3(floor(size(model(von3:bis3,3),1)/2)+1-ftc:floor(size(model(von3:bis3,3),1)/2)+ftc) = zeros(ftc*2,1);
model(von3:bis3,:) = [real(ifft(f1)) real(ifft(f2)) real(ifft(f3))];


figure
hold off
plot3(model(:,1),model(:,2),model(:,3),'r')
axis equal

if ~(size(unten_in,1) == size(unten,1) && size(unten,1) == size(unten_legs,1))
    disp('Problem! Size of all arrays should be the same!!')
end

x = [ repmat(model(1,1),1,size(unten,1)); ...
        [model(von3:bis3,1)',model(von3,1)']; ...
        [model(von2:bis2,1)',model(von2,1)']; ...
        [model(von1:bis1,1)',model(von1,1)']; ...
        [model(von:bis,1)',model(von,1)']; ...
        [model(von:bis,1)',model(von,1)'] ];
y = [ repmat(model(1,2),1,size(unten,1)); ...
        [model(von3:bis3,2)',model(von3,2)']; ...
        [model(von2:bis2,2)',model(von2,2)']; ...
        [model(von1:bis1,2)',model(von1,2)']; ...
        [model(von:bis,2)',model(von,2)']; ...
        [model(von:bis,2)',model(von,2)'] ];
z = [ repmat(model(1,3),1,size(unten,1)); ...
        [model(von3:bis3,3)',model(von3,3)']; ...
        [model(von2:bis2,3)',model(von2,3)']; ...
        [model(von1:bis1,3)',model(von1,3)']; ...
        [model(von:bis,3)',model(von,3)']; ...
        zeros(1,size(unten,1)) ];

hold on
surf(x,y,z)
axis equal
colormap(lines)
drawnow

disp('writing chamber model1 ...')
surf2stl(c_file,x,y,z)


%% chamber model 2

figure

model2 = [1 0 0; 0 cos(chamber_angle) -sin(chamber_angle); 0 sin(chamber_angle) cos(chamber_angle)] ...
        * oben(1,:)';

model2 = [model2 [1 0 0; 0 cos(chamber_angle) -sin(chamber_angle); 0 sin(chamber_angle) cos(chamber_angle)] ...
        * oben_in(2:size(oben_in,1),:)' - repmat(model2(:,1),1,size(oben_in,1)-1) ];
model2 = [model2 [1 0 0; 0 cos(chamber_angle) -sin(chamber_angle); 0 sin(chamber_angle) cos(chamber_angle)] ...
        * oben(2:size(oben,1),:)' - repmat(model2(:,1),1,size(oben,1)-1) ];
model2 = [model2 [1 0 0; 0 cos(chamber_angle) -sin(chamber_angle); 0 sin(chamber_angle) cos(chamber_angle)] ...
        * unten(2:size(unten,1),:)' - repmat(model2(:,1)-[0;0;chamber_in_skull+chamber_leg_thickness],1,size(unten,1)-1) ];
model2 = [model2 [1 0 0; 0 cos(chamber_angle) -sin(chamber_angle); 0 sin(chamber_angle) cos(chamber_angle)] ...
        * unten_legs(2:size(unten_legs,1),:)' - repmat(model2(:,1)-[0;0;chamber_leg_thickness],1,size(unten_legs,1)-1) ];
model2(:,1) = [0;0;-3];
    
model2 = model2';

%fft smooth unten

bis0 = size(model2,1);
von0 = 1+bis0 - (size(unten_legs,1)-1);
f1 = fft(model2(von0:bis0,1));
f2 = fft(model2(von0:bis0,2));
f3 = fft(model2(von0:bis0,3));
f1(floor(size(model2(von0:bis0,1),1)/2)+1-ftc:floor(size(model2(von0:bis0,1),1)/2)+ftc) = zeros(ftc*2,1);
f2(floor(size(model2(von0:bis0,2),1)/2)+1-ftc:floor(size(model2(von0:bis0,2),1)/2)+ftc) = zeros(ftc*2,1);
f3(floor(size(model2(von0:bis0,3),1)/2)+1-ftc:floor(size(model2(von0:bis0,3),1)/2)+ftc) = zeros(ftc*2,1);
model2(von0:bis0,:) = [real(ifft(f1)) real(ifft(f2)) real(ifft(f3))];

bis01 = von0-1;
von01 = 1+bis01 - (size(unten,1)-1);
f1 = fft(model2(von01:bis01,1));
f2 = fft(model2(von01:bis01,2));
f3 = fft(model2(von01:bis01,3));
f1(floor(size(model2(von01:bis01,1),1)/2)+1-ftc:floor(size(model2(von01:bis01,1),1)/2)+ftc) = zeros(ftc*2,1);
f2(floor(size(model2(von01:bis01,2),1)/2)+1-ftc:floor(size(model2(von01:bis01,2),1)/2)+ftc) = zeros(ftc*2,1);
f3(floor(size(model2(von01:bis01,3),1)/2)+1-ftc:floor(size(model2(von01:bis01,3),1)/2)+ftc) = zeros(ftc*2,1);
model2(von01:bis01,:) = [real(ifft(f1)) real(ifft(f2)) real(ifft(f3))];

bis02 = von01-1;
von02 = 1+bis02 - (size(oben,1)-1);

bis03 = von02-1;
von03 = 1+bis03 - (size(oben_in,1)-1);

model2(:,3) = model2(:,3) - min(model2(:,3)) + chamber_leg_thickness;


hold off
plot3(model2(:,1),model2(:,2),model2(:,3),'r')
axis equal

if ~(size(oben,1) == size(oben_in,1) && size(unten,1) == size(unten_legs,1))
    disp('Problem! Size of all arrays should be the same!!')
end

x = [ repmat(model2(1,1),1,size(oben,1)); ...
        [model2(von03:bis03,1)',model2(von03,1)']; ...
        [model2(von02:bis02,1)',model2(von02,1)']; ...
        [model2(von01:bis01,1)',model2(von01,1)']; ...
        [model2(von0:bis0,1)',model2(von0,1)']; ...
        [model2(von0:bis0,1)',model2(von0,1)'] ];
y = [ repmat(model2(1,2),1,size(oben,1)); ...
        [model2(von03:bis03,2)',model2(von03,2)']; ...
        [model2(von02:bis02,2)',model2(von02,2)']; ...
        [model2(von01:bis01,2)',model2(von01,2)']; ...
        [model2(von0:bis0,2)',model2(von0,2)']; ...
        [model2(von0:bis0,2)',model2(von0,2)'] ];
z = [ repmat(model2(1,3),1,size(oben,1)); ...
        [model2(von03:bis03,3)',model2(von03,3)']; ...
        [model2(von02:bis02,3)',model2(von02,3)']; ...
        [model2(von01:bis01,3)',model2(von01,3)']; ...
        [model2(von0:bis0,3)',model2(von0,3)']; ...
        zeros(1,size(unten_legs,1)) ];

hold on
surf(x,y,z)
axis equal
colormap(lines)
drawnow

disp('writing chamber model2 ...')
surf2stl(c_file2,x,y,z)

end

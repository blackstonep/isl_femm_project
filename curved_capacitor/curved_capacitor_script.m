clear
% Analyze capicitance of two parallel plates for
% configuration 1: lipped bottom plate and
% 2: sphereical or curved bottom plate.
mu_0 = 4*pi*1e-7;
eps_0 = 8.8542e-15; % units of millimeters. 
top_plate_voltage = 1.0;



dmax = 0.015;    % d = separation of plates
R = 40;            % radius of plates
hmax = 0.014;    % h = depth difference of bottom plate
hmin = -0.2;
thick = 0.040;    % thickness of plates
Number_of_Steps = 1000;

% Prospect for loop:
d = dmax;
h = hmax;

%Positions of top plate:
top_llx = 0.0;
top_lly = d;
top_ulx = 0.0;
top_uly = d + thick;
top_lrx = R;
top_lry = d;
top_urx = R;
top_ury = d + thick;

%...and bottom plate:
bot_urx = R;
bot_ury = h;
bot_lrx = R;
bot_lry = h-thick;
bot_ulx = 0.0;
bot_uly = 0.0;
bot_llx = 0.0;
bot_lly = 0.0-thick;

%Bottom plate joint nodes:
h_off = 0.1;                % x-dir offset for slope
bot_joint_lx = R - h_off;
bot_joint_ly = 0.0 - thick;
bot_joint_ux = R - h_off;
bot_joint_uy = 0.0;

%Positions of Vacuum can. 
can_llx = 0;
can_lly = -110;
can_lrx = 90;
can_lry = -110;
can_urx = 90;
can_ury = 110;
can_ulx = 0;
can_uly = 110;

%Calculate arcmeasure of lower plate in degrees:
    %Ah! Femm only takes 1 < arc < 180. :(
%alpha = asin(R / sqrt(R*R + h*h));
%beta = asin(h / sqrt(R*R + h*h));
%arcm = (pi/2.0 - alpha + beta)*180.0 / pi;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

openfemm

opendocument('polygon_capacitor.FEE');
ei_probdef('millimeters','axi',1e-8,0,30)
% automatic mesh sizes [call to triangle?]
    air_mesh = 1;
    cu_mesh = 1;

%Draw top capacitor plate:
ei_addnode(top_llx, top_lly);
ei_addnode(top_ulx, top_uly);
ei_addnode(top_lrx, top_lry);
ei_addnode(top_urx, top_ury);
ei_addsegment(top_llx, top_lly, top_lrx, top_lry);
ei_addsegment(top_lrx, top_lry, top_urx, top_ury);
ei_addsegment(top_urx, top_ury, top_ulx, top_uly);
ei_clearselected;

%Draw bottom capacitor plate:
ei_addnode(bot_llx, bot_lly);
ei_addnode(bot_ulx, bot_uly);
ei_addnode(bot_lrx, bot_lry);
ei_addnode(bot_urx, bot_ury);
ei_addnode(bot_joint_lx, bot_joint_ly);
ei_addnode(bot_joint_ux, bot_joint_uy);
ei_addsegment(bot_llx, bot_lly, bot_ulx, bot_uly);
ei_addsegment(bot_ulx, bot_uly, bot_joint_ux, bot_joint_uy);
ei_addsegment(bot_joint_ux, bot_joint_uy, bot_urx, bot_ury);
ei_addsegment(bot_urx, bot_ury, bot_lrx, bot_lry);
ei_addsegment(bot_lrx, bot_lry, bot_joint_lx, bot_joint_ly);
ei_addsegment(bot_joint_lx, bot_joint_ly, bot_llx, bot_lly);
ei_clearselected;

%Draw vacuum can.
ei_addnode(can_llx, can_lly);
ei_addnode(can_lrx, can_lry);
ei_addnode(can_urx, can_ury);
ei_addnode(can_ulx, can_uly);
ei_addsegment(can_llx, can_lly, can_lrx, can_lry);
ei_addsegment(can_lrx, can_lry, can_urx, can_ury);
ei_addsegment(can_urx, can_ury, can_ulx, can_uly);
ei_addsegment(can_ulx, can_uly, can_llx, can_lly);
ei_clearselected;

% Coordinate labels:
    %air
    x_air = R / 2.0;
    y_air = R / 4.0;

    %Upper plate
    x_upper = R / 2.0;
    y_upper = d + 0.5*thick;

    %Lower plate
    x_lower = R / 2.0;
    y_lower = 0.0-0.5*thick;

%Define material properties
ei_addmaterial('Air', 1, 1, 0);
ei_addmaterial('Copper', 1, 1, 0);

% Label air as group 3
ei_addblocklabel(x_air, y_air);
ei_selectlabel(x_air, y_air);
ei_setblockprop('Air', air_mesh, 2, 3);
ei_clearselected;

%Top plate as group 1
ei_addblocklabel(x_upper, y_upper);
ei_selectlabel(x_upper, y_upper);
ei_setblockprop('Copper', cu_mesh, 2, 1);
ei_clearselected;

%Lower plate as group 2
ei_addblocklabel(x_lower, y_lower);
ei_selectlabel(x_lower, y_lower);
ei_setblockprop('Copper', cu_mesh, 2, 2);
ei_clearselected;


% Create boundary property for the capacitors
ei_addboundprop('top_plate_boundary', top_plate_voltage, 0, 0, 0, 0);
ei_addboundprop('bot_plate_boundary',0.0, 0, 0, 0, 0);

% Apply boundary properties

% Top Plate
ei_selectsegment(0, 0.5*(top_lly+top_uly));
ei_selectsegment(R/2.0, d);
ei_selectsegment(R, 0.5*(top_lry+top_ury));
ei_selectsegment(R/2.0, d + thick);
ei_setsegmentprop('top_plate_boundary', 0.1, 1, 0, 1,'<None>');
ei_clearselected;
    %Add nodes to group as well for cleaning!
    ei_selectnode(0, d);
    ei_selectnode(R, d);
    ei_selectnode(R, d + thick);
    ei_selectnode(0, d + thick);
    ei_setnodeprop('<None>', 1, '<None>');
    ei_clearselected;

% Bottom Plate
ei_selectsegment(0, -0.5*thick);
ei_selectsegment(R/2, 0);
ei_selectsegment(R/2, -thick);
ei_selectsegment(R - 0.5*h_off, 0.5*h);
ei_selectsegment(R - 0.5*h_off, 0.5*h - thick);
ei_selectsegment(R, 0.5*thick - h);
ei_setsegmentprop('bot_plate_boundary', 0.1, 1, 0, 2, '<None>');
ei_clearselected;
    %Add nodes to group as well for cleaning!
    ei_selectnode(0,0);
    ei_selectnode(0, -thick);
    ei_selectnode(R, h);
    ei_selectnode(R - h_off, 0);
    ei_selectnode(R, h - thick);
    ei_selectnode(R - h_off, -thick);
    ei_setnodeprop('<None>', 2, '<None>');
    ei_clearselected;

%Vacuum can, outside edges grounded to bottom plate. 
ei_selectsegment(R/2, can_lly);
ei_selectsegment(can_lrx, 0);
ei_selectsegment(R/2, can_uly);
ei_setsegmentprop('bot_plate_boundary',0.1,1,0,3,'<None>');
ei_clearselected;



nd = Number_of_Steps; %Number of iterations in loop:
hinc = (hmax - hmin) / nd; %Incrementation per iterations

format shortEng

for j = 0:nd

    % Runs the electrostatic solver:
    ei_analyze();

    ei_loadsolution;

    % Record energy stored in capacitor
    eo_selectblock(R/2, R/4);
    Integral = eo_blockintegral(0);
    Energy_num = Integral(1);
    Energy_theory = 0.5*eps_0*pi*R*R / d * (top_plate_voltage)*(top_plate_voltage);
    Relative_Error = abs(Energy_num - Energy_theory)/Energy_theory;

    h_new = bot_ury - j*hinc;

    zeta = bot_joint_lx + R*d/h_new;
    log_arg = 1 + bot_joint_lx*h_new/(R*d) - h_new/d;
    theory_capacitance = 2*pi*eps_0*(bot_joint_lx^2/(2*d) - R/h_new * (R - bot_joint_lx + zeta*log(log_arg)) );

    ana_nrg(j+1) = 0.5*theory_capacitance*(top_plate_voltage^2);
    s(j+1) = Energy_num;
    H(j+1) = hmax - hinc*j;
    err(j+1) = Relative_Error;
    ana_nrg_err(j+1) = abs(ana_nrg(j+1)-s(j+1))/ana_nrg(j+1);

    eo_clearblock;
    
    ei_selectnode(R, hmax - hinc*j); 
    ei_selectnode(R, hmax - hinc*j - thick);    
    ei_movetranslate(0, -hinc);
    ei_clearselected;   


    data_matrix(1,j+1) = H(j+1);
    data_matrix(2,j+1) = s(j+1);
    data_matrix(3,j+1) = err(j+1);
    data_matrix(4,j+1) = ana_nrg(j+1);
    data_matrix(5,j+1) = ana_nrg_err(j+1);

end 

fid = fopen('poly.dat', 'w');
fprintf(fid, '%50s\n', '#===============================================================================');
fprintf(fid, '%s %12s %15s %15s %15s %15s\n', '#', 'Lip Height', 'Numerical Energy', 'Relative Error', 'Analytical Energy', 'Analyt. Enrgy. Error');
fprintf(fid, '%50s\n', '#===============================================================================');
fprintf(fid, '%12.8f %.15f %12.8f %.15f %12.8f\n', data_matrix);


%Plot stored energy versus lip size.
%figure(1)
%plot(H, s)
%xlabel('Height of Lip above plate capacitor (mm)')
%ylabel('Stored Potential Energy (J)')

%Plot deviation with height of lip.
%figure(2)
%plot(H, err)
%xlabel('Height of Lip above plate capacitor (mm)')
%ylabel('Relatvie Deviation of Num. Energy to Expected (J)')
%title('Relative Deviation from Ideal, flat Plate Capacitors')

%Erase previous plates
ei_selectgroup(1);
ei_deleteselectednodes;
ei_deleteselected;
ei_selectgroup(2);
ei_deleteselectednodes;
ei_deleteselected;

status = fclose(fid);

%ei_saveas();
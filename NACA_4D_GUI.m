function varargout = NACA_4D_GUI(varargin)
% NACA_4D_GUI MATLAB code for NACA_4D_GUI.fig
%      NACA_4D_GUI, by itself, creates a new NACA_4D_GUI or raises the existing
%      singleton*.
%
%      H = NACA_4D_GUI returns the handle to a new NACA_4D_GUI or the handle to
%      the existing singleton*.
%
%      NACA_4D_GUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in NACA_4D_GUI.M with the given input arguments.
%
%      NACA_4D_GUI('Property','Value',...) creates a new NACA_4D_GUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before NACA_4D_GUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to NACA_4D_GUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help NACA_4D_GUI

% Last Modified by GUIDE v2.5 05-Dec-2017 08:32:34

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @NACA_4D_GUI_OpeningFcn, ...
                   'gui_OutputFcn',  @NACA_4D_GUI_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before NACA_4D_GUI is made visible.
function NACA_4D_GUI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to NACA_4D_GUI (see VARARGIN)

% Choose default command line output for NACA_4D_GUI
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes NACA_4D_GUI wait for user response (see UIRESUME)
% uiwait(handles.figure1);

% --- Outputs from this function are returned to the command line.
function varargout = NACA_4D_GUI_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


function Airfoil11_Callback(hObject, eventdata, handles)
% hObject    handle to Airfoil11 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Airfoil11 as text
% str2double(get(hObject,'String')) returns contents of Airfoil11 as a double


% --- Executes during object creation, after setting all properties.
function Airfoil11_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Airfoil11 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
% See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function Alpha1_Callback(hObject, eventdata, handles)
% hObject    handle to Alpha1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Alpha1 as text
% str2double(get(hObject,'String')) returns contents of Alpha1 as a double


% --- Executes during object creation, after setting all properties.
function Alpha1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Alpha1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
% See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in Chord1.
function Chord1_Callback(hObject, eventdata, handles)
% hObject    handle to Chord1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of Chord1


% --- Executes on button press in MeanCamberLine1.
function MeanCamberLine1_Callback(hObject, eventdata, handles)
% hObject    handle to MeanCamberLine1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of MeanCamberLine1


% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
clc
MCL = get(handles.MeanCamberLine1,'Value');
radial_density = str2double(get(handles.Rad_density ,'string'));
gridpoints = str2double(get(handles.Num_points ,'string'));

Airfoil = get(handles.Airfoil11 ,'string');
chordlength = 1;
%gridpoints = 300;
%radial_density = 50;
AngleofAttack = str2double(get(handles.Alpha1,'string'));  %degrees

%% Airfoil Drawing
disp('Drawing Airfoil')

m  = str2double(Airfoil(1));
p  = str2double(Airfoil(2));
xx = str2double(Airfoil(3:4));

%define parameters
M = m/100;
P = p/10;
T = xx/100;
c = chordlength;
step = c/gridpoints;

%Mean Camber Line:
x = zeros(gridpoints,1);
for i = 1:gridpoints
    x(i) = 0.5*cos(2*pi*i/gridpoints)+0.5;
end
x = sort(x);

%x = [0:step:c].';
y_c = [0];
dy_c = [0];
for j = 1:gridpoints
   if x(j) >=0 && x(j) < P.*c;
        y_c(j,1)  = M/(P.^2).*x(j).*(2.*P - x(j)/c);
        dy_c(j,1) = 2.*M/(P.^2).*(P - x(j)/c);     
   else x(j) >= P.*c && x(j) <= c;
        y_c(j,1)  = M.*((c-x(j))/((1-P).^2)).*( 1+x(j)/c - 2.*P );
        dy_c(j,1) = 2.*M/((1 - P).^2).*(x(j)/c - P);  
   end
end

%Angle of Attack 
alpha = -AngleofAttack*pi/180;
Rotation = ([cos(alpha), -sin(alpha); sin(alpha), cos(alpha)]*[x,y_c].').';
x = Rotation(:,1);
y_c = Rotation(:,2);

%Thinkness Distribution
a_0 = 0.2969;
a_1 = -0.126;
a_2 = -0.3516;
a_3 = 0.2843;
a_4 = -0.1015;   %trailing edge pointy
%a_4 = -0.1036;   %trailing edge closed

y_t = T/0.2.*(  a_0.*x.^(0.5) + a_1.*x + a_2.*x.^2 + a_3.*x.^3 + a_4.*x.^4);

theta = atan(dy_c);
x_u = x - y_t.*sin(theta)-0.1;
y_u = y_c + y_t.*cos(theta);

x_L = x + y_t.*sin(theta)-0.1;
y_L = y_c - y_t.*cos(theta);
hold off
plot(0,0)
hold on
axis equal
plot(x_u,y_u,'r-')
plot(x_L,y_L,'r-')

if MCL == 1
plot(x-0.1,y_c,'b-')
else
    %Do nothing
end

function Num_points_Callback(hObject, eventdata, handles)
% hObject    handle to Num_points (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Num_points as text
% str2double(get(hObject,'String')) returns contents of Num_points as a double


% --- Executes during object creation, after setting all properties.
function Num_points_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Num_points (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
% See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function Rad_density_Callback(hObject, eventdata, handles)
% hObject    handle to Rad_density (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Rad_density as text
% str2double(get(hObject,'String')) returns contents of Rad_density as a double


% --- Executes during object creation, after setting all properties.
function Rad_density_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Rad_density (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
% See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in Gen_Mesh.
function Gen_Mesh_Callback(hObject, eventdata, handles)
% hObject    handle to Gen_Mesh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

hold off
plot(0,0)

MCL = get(handles.MeanCamberLine1,'Value');
radial_density = str2double(get(handles.Rad_density ,'string'));
gridpoints = str2double(get(handles.Num_points ,'string'));

Airfoil = get(handles.Airfoil11 ,'string');
chordlength = 1;
%gridpoints = 300;
%radial_density = 50;
AngleofAttack = str2double(get(handles.Alpha1,'string'));  %degrees

%% Airfoil Drawing
disp('Drawing Airfoil')

m  = str2double(Airfoil(1));
p  = str2double(Airfoil(2));
xx = str2double(Airfoil(3:4));

%define parameters
M = m/100;
P = p/10;
T = xx/100;
c = chordlength;
step = c/gridpoints;


%Mean Camber Line:
x = zeros(gridpoints,1);
for i = 1:gridpoints
    x(i) = 0.5*cos(2*pi*i/gridpoints)+0.5;
end
x = sort(x);

%x = [0:step:c].';
y_c = [0];
dy_c = [0];
for j = 1:gridpoints
   if x(j) >=0 && x(j) < P.*c;
        y_c(j,1)  = M/(P.^2).*x(j).*(2.*P - x(j)/c);
        dy_c(j,1) = 2.*M/(P.^2).*(P - x(j)/c);    
   else x(j) >= P.*c && x(j) <= c;
        y_c(j,1)  = M.*((c-x(j))/((1-P).^2)).*( 1+x(j)/c - 2.*P );
        dy_c(j,1) = 2.*M/((1 - P).^2).*(x(j)/c - P);   
   end
end


%Angle of Attack 
alpha = -AngleofAttack*pi/180;
Rotation = ([cos(alpha), -sin(alpha); sin(alpha), cos(alpha)]*[x,y_c].').';
x = Rotation(:,1);
y_c = Rotation(:,2);

hold on

%Thinkness Distribution
a_0 = 0.2969;
a_1 = -0.126;
a_2 = -0.3516;
a_3 = 0.2843;
a_4 = -0.1015;   %trailing edge pointy
%a_4 = -0.1036;  %trailing edge closed

y_t = T/0.2.*(  a_0.*x.^(0.5) + a_1.*x + a_2.*x.^2 + a_3.*x.^3 + a_4.*x.^4);

theta = atan(dy_c);
x_u = x - y_t.*sin(theta)-0.1;
y_u = y_c + y_t.*cos(theta);

x_L = x + y_t.*sin(theta)-0.1;
y_L = y_c - y_t.*cos(theta);

%% Mesh Point Generation 
disp('Generating Mesh')

%Start the point matrix with the airfoil points
Points(1:gridpoints,1) = x_u;
Points(1:gridpoints,2) = y_u;
Points( gridpoints+1:2*(gridpoints) , 1 ) = flipud(x_L);
Points( gridpoints+1:2*(gridpoints) , 2 ) = flipud(y_L);
enis = length(Points);

for i = 1:enis
        for j = 1:radial_density
            Points(2*(gridpoints) + radial_density*i-(radial_density-j),:) = (exp(2.8*j/radial_density)).*Points(i,:); 
        end
end

%Throw away unnessisary points
sqweenis = length(Points);
for i = 1:sqweenis
    if abs(Points(sqweenis-i+1,1)-0.5)> 1.05
        Points(sqweenis-i+1,:) = [];
    elseif abs( Points(sqweenis-i+1,2) ) > 0.55
        Points(sqweenis-i+1,:) = [];
    end
end

DT = delaunayTriangulation(Points);
triplot(DT,'k')

if MCL == 1
plot(x-0.1,y_c,'b-')
else
    %Do nothing
end

plot(x_u,y_u,'r-')
plot(x_L,y_L,'r-')
axis([-0.5,1.5,-0.5,0.5])
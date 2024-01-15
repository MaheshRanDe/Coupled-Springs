function varargout = Coupled_springs(varargin)
% COUPLED_SPRINGS MATLAB code for Coupled_springs.fig
%      COUPLED_SPRINGS, by itself, creates a new COUPLED_SPRINGS or raises the existing
%      singleton*.
%
%      H = COUPLED_SPRINGS returns the handle to a new COUPLED_SPRINGS or the handle to
%      the existing singleton*.
%
%      COUPLED_SPRINGS('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in COUPLED_SPRINGS.M with the given input arguments.
%
%      COUPLED_SPRINGS('Property','Value',...) creates a new COUPLED_SPRINGS or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before Coupled_springs_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to Coupled_springs_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help Coupled_springs

% Last Modified by GUIDE v2.5 25-Apr-2019 03:58:15

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @Coupled_springs_OpeningFcn, ...
                   'gui_OutputFcn',  @Coupled_springs_OutputFcn, ...
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


% --- Executes just before Coupled_springs is made visible.
function Coupled_springs_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to Coupled_springs (see VARARGIN)

% Choose default command line output for Coupled_springs
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes Coupled_springs wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = Coupled_springs_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



function edit1_Callback(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit1 as text
%        str2double(get(hObject,'String')) returns contents of edit1 as a double


% --- Executes during object creation, after setting all properties.
function edit1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit2_Callback(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit2 as text
%        str2double(get(hObject,'String')) returns contents of edit2 as a double


% --- Executes during object creation, after setting all properties.
function edit2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit3_Callback(hObject, eventdata, handles)
% hObject    handle to edit3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit3 as text
%        str2double(get(hObject,'String')) returns contents of edit3 as a double


% --- Executes during object creation, after setting all properties.
function edit3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit4_Callback(hObject, eventdata, handles)
% hObject    handle to edit4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit4 as text
%        str2double(get(hObject,'String')) returns contents of edit4 as a double


% --- Executes during object creation, after setting all properties.
function edit4_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit5_Callback(hObject, eventdata, handles)
% hObject    handle to edit5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit5 as text
%        str2double(get(hObject,'String')) returns contents of edit5 as a double


% --- Executes during object creation, after setting all properties.
function edit5_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit6_Callback(hObject, eventdata, handles)
% hObject    handle to edit6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit6 as text
%        str2double(get(hObject,'String')) returns contents of edit6 as a double


% --- Executes during object creation, after setting all properties.
function edit6_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit7_Callback(hObject, eventdata, handles)
% hObject    handle to edit7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit7 as text
%        str2double(get(hObject,'String')) returns contents of edit7 as a double


% --- Executes during object creation, after setting all properties.
function edit7_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit8_Callback(hObject, eventdata, handles)
% hObject    handle to edit8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit8 as text
%        str2double(get(hObject,'String')) returns contents of edit8 as a double


% --- Executes during object creation, after setting all properties.
function edit8_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit9_Callback(hObject, eventdata, handles)
% hObject    handle to edit9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit9 as text
%        str2double(get(hObject,'String')) returns contents of edit9 as a double


% --- Executes during object creation, after setting all properties.
function edit9_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit10_Callback(hObject, eventdata, handles)
% hObject    handle to edit10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit10 as text
%        str2double(get(hObject,'String')) returns contents of edit10 as a double


% --- Executes during object creation, after setting all properties.
function edit10_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit11_Callback(hObject, eventdata, handles)
% hObject    handle to edit11 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit11 as text
%        str2double(get(hObject,'String')) returns contents of edit11 as a double


% --- Executes during object creation, after setting all properties.
function edit11_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit11 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit12_Callback(hObject, eventdata, handles)
% hObject    handle to edit12 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit12 as text
%        str2double(get(hObject,'String')) returns contents of edit12 as a double


% --- Executes during object creation, after setting all properties.
function edit12_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit12 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit13_Callback(hObject, eventdata, handles)
% hObject    handle to edit13 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit13 as text
%        str2double(get(hObject,'String')) returns contents of edit13 as a double


% --- Executes during object creation, after setting all properties.
function edit13_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit13 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit14_Callback(hObject, eventdata, handles)
% hObject    handle to edit14 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit14 as text
%        str2double(get(hObject,'String')) returns contents of edit14 as a double


% --- Executes during object creation, after setting all properties.
function edit14_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit14 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit15_Callback(hObject, eventdata, handles)
% hObject    handle to edit15 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit15 as text
%        str2double(get(hObject,'String')) returns contents of edit15 as a double


% --- Executes during object creation, after setting all properties.
function edit15_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit15 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit16_Callback(hObject, eventdata, handles)
% hObject    handle to edit16 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit16 as text
%        str2double(get(hObject,'String')) returns contents of edit16 as a double


% --- Executes during object creation, after setting all properties.
function edit16_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit16 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit17_Callback(hObject, eventdata, handles)
% hObject    handle to edit17 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit17 as text
%        str2double(get(hObject,'String')) returns contents of edit17 as a double


% --- Executes during object creation, after setting all properties.
function edit17_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit17 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit18_Callback(hObject, eventdata, handles)
% hObject    handle to edit18 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit18 as text
%        str2double(get(hObject,'String')) returns contents of edit18 as a double


% --- Executes during object creation, after setting all properties.
function edit18_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit18 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit19_Callback(hObject, eventdata, handles)
% hObject    handle to edit19 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit19 as text
%        str2double(get(hObject,'String')) returns contents of edit19 as a double


% --- Executes during object creation, after setting all properties.
function edit19_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit19 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit20_Callback(hObject, eventdata, handles)
% hObject    handle to edit20 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit20 as text
%        str2double(get(hObject,'String')) returns contents of edit20 as a double


% --- Executes during object creation, after setting all properties.
function edit20_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit20 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%extracting values of paremeter
m1=str2num(get(handles.edit1,'String'));
m2=str2num(get(handles.edit2,'String'));
L1=str2num(get(handles.edit3,'String'));
L2=str2num(get(handles.edit4,'String'));
k1=str2num(get(handles.edit5,'String'));
k2=str2num(get(handles.edit6,'String'));
x1=str2num(get(handles.edit13,'String'));
x2=str2num(get(handles.edit14,'String'));
v1=str2num(get(handles.edit15,'String'));
v2=str2num(get(handles.edit16,'String'));
T=str2num(get(handles.edit22,'String'));
h=str2num(get(handles.edit21,'String'));
g=str2num(get(handles.edit23,'String'));
I=str2num(get(handles.radiobutton1,'String'));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
y0=[x1;x2;v1;v2];%initial value
r=0.01;%radius for spring
t0=0;n=(T-t0)/h;t=linspace(t0,T,n);%assignining time values and calculating number os steps to be iterated
f=@(t,x) [x(3);x(4);-k1/m1*x(1)-k2/m1*(x(1)-x(2));-k2/m2*(x(2)-x(1))];%defining the function
x=zeros(4,n);x(:,1)=y0;%defining the matrix for aproximations and initializing
[t1,Y]=ode45(f,t,y0);%calculating exact values
%defining RK-method
A=[0 0 0;0.5 0 0;-1 2 0];
c=[0 0.5 1];
b=[1/6 2/3 1/6];
xi=zeros(4,3);%RK-stages
xi(:,1)=y0;
for i=1:n-1
    if I==0
    xi(:,1)=y0;
    xi(:,2)=y0+h*A(2,1)*f(t(i)+c(1)*h,xi(:,1));
    xi(:,3)=y0+h*A(3,1)*f(t(i)+c(1)*h,xi(:,1))+h*A(3,2)*f(t(i)+c(2)*h,xi(:,2));
    y0=y0+h*(b(1)*f(t(i)+c(1)*h,xi(:,1))+b(2)*f(t(i)+c(2)*h,xi(:,2))+b(3)*f(t(i)+c(3)*h,xi(:,3)));
    x(:,i+1)=y0;
    else
        p=@(y) y-x(:,i)-h*f((t(i+1)+t(i))/2,(x(:,i)+y)/2);
       x(:,i+1)=fsolve(@(y)p(y),x(:,i));
    end
end
axes(handles.axes2)
 plot(t,x(1,:),'r-',t,Y(:,1),'b-');
 xlabel('Time(s)');
 ylabel('Elongation of spring 1');
 legend('Aproximation','Exact');
 title('Elongation of spring 1 vs time')
 axes(handles.axes3)
 plot(t,x(2,:),'r-',t,Y(:,2),'b-');
 xlabel('Time(s)');
 ylabel('Elongation of spring 2');
 legend('Aproximation','Exact');
 title('Elongation of spring 2 vs time');
for i=1:n-1
    axes(handles.axes1)
%%
tt = pi/2:pi/50:21*pi/2;
xx=linspace(0,L1+x(1,i),length(tt));
xx2=linspace(L1+x(1,i),L1+L2+x(1,i)+x(2,i),length(tt));
[Xs,Ys,Zs]=sphere(50);
st =0.01*sin(tt);
ct =0.01*cos(tt);
xmax=max(L1+L2+x(1,:)+x(2,:))+0.2;
%%
xplane=linspace(0,xmax,10);yplane=linspace(-0.1,0.1,10);
[Xplane,Yplane]=meshgrid(xplane,yplane);
[Yplane1,Zplane1]=meshgrid(yplane,yplane);
Zplane=-0.01*ones(size(Xplane));
Xplane1=zeros(size(Yplane));
Yplane2=0.1*ones(size(Yplane));
surf(Xplane,Yplane,Zplane,'FaceColor',[0 0.5 0],'EdgeColor','None');
hold on
surf(Xplane1,Yplane1,Zplane1,'FaceColor',[0 1 1],'EdgeColor','None');
surf(Xplane,Yplane2,Zplane1,'FaceColor',[0 1 1],'EdgeColor','None');
view(60, 45)
light
lighting phong
%%
plot3(xx,st,ct,'r')
plot3(xx2,st,ct,'r')
%%
surf(L1+x(1,i)+r*Xs,r*Ys,r*Zs,'FaceColor',[0 0 1],'EdgeColor','None');
surf(L1+L2+x(1,i)+x(2,i)+r*Xs,r*Ys,r*Zs,'FaceColor',[0 0 1],'EdgeColor','None');
axis([0 xmax -0.1 0.1 -0.01 0.1])
xlabel('x');ylabel('y');zlabel('z');
getframe;
hold off
end


% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%extracting values of paremeter
m1=str2num(get(handles.edit1,'String'));
m2=str2num(get(handles.edit2,'String'));
L1=str2num(get(handles.edit3,'String'));
L2=str2num(get(handles.edit4,'String'));
k1=str2num(get(handles.edit5,'String'));
k2=str2num(get(handles.edit6,'String'));
x1=str2num(get(handles.edit13,'String'));
x2=str2num(get(handles.edit14,'String'));
v1=str2num(get(handles.edit15,'String'));
v2=str2num(get(handles.edit16,'String'));
T=str2num(get(handles.edit22,'String'));
h=str2num(get(handles.edit21,'String'));
g=str2num(get(handles.edit23,'String'));
mu1=str2num(get(handles.edit7,'String'));
mu2=str2num(get(handles.edit8,'String'));
alpha1=str2num(get(handles.edit9,'String'));
alpha2=str2num(get(handles.edit10,'String'));
beta1=str2num(get(handles.edit11,'String'));
beta2=str2num(get(handles.edit12,'String'));
I=str2num(get(handles.radiobutton1,'String'));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
y0=[x1;x2;v1;v2];%initial value
r=0.01;%radius for spring
t0=0;n=(T-t0)/h;t=linspace(t0,T,n);%assignining time values and calculating number os steps to be iterated
f=@(t,x) [x(3);x(4);-mu1*g-alpha1*x(3)+beta1*x(1)^3+beta2*(x(1)-x(2))^3-k1/m1*x(1)-k2/m1*(x(1)-x(2));-mu2*g-alpha2*x(4)+beta2*(x(2)-x(1))^3-k2/m2*(x(2)-x(1))];%defining the function
x=zeros(4,n);x(:,1)=y0;%defining the matrix for aproximations and initializing
[t1,Y]=ode45(f,t,y0);%calculating exact values
%defining RK-method
A=[0 0 0;0.5 0 0;-1 2 0];
c=[0 0.5 1];
b=[1/6 2/3 1/6];
xi=zeros(4,3);%RK-stages
xi(:,1)=y0;
for i=1:n-1
    if I==0
    xi(:,1)=y0;
    xi(:,2)=y0+h*A(2,1)*f(t(i)+c(1)*h,xi(:,1));
    xi(:,3)=y0+h*A(3,1)*f(t(i)+c(1)*h,xi(:,1))+h*A(3,2)*f(t(i)+c(2)*h,xi(:,2));
    y0=y0+h*(b(1)*f(t(i)+c(1)*h,xi(:,1))+b(2)*f(t(i)+c(2)*h,xi(:,2))+b(3)*f(t(i)+c(3)*h,xi(:,3)));
    x(:,i+1)=y0;
    else
        p=@(y) y-x(:,i)-h*f((t(i+1)+t(i))/2,(x(:,i)+y)/2);
       x(:,i+1)=fsolve(@(y)p(y),x(:,i));
    end
end
axes(handles.axes2)
 plot(t,x(1,:),'r-',t,Y(:,1),'b-');
 xlabel('Time(s)');
 ylabel('Elongation of spring 1');
 legend('Aproximation','Exact');
 title('Elongation of spring 1 vs time')
 axes(handles.axes3)
 plot(t,x(2,:),'r-',t,Y(:,2),'b-');
 xlabel('Time(s)');
 ylabel('Elongation of spring 2');
 legend('Aproximation','Exact');
 title('Elongation of spring 2 vs time');
for i=1:n-1
    axes(handles.axes1)
%%
tt = pi/2:pi/50:21*pi/2;
xx=linspace(0,L1+x(1,i),length(tt));
xx2=linspace(L1+x(1,i),L1+L2+x(1,i)+x(2,i),length(tt));
[Xs,Ys,Zs]=sphere(50);
st =0.01*sin(tt);
ct =0.01*cos(tt);
xmax=max(L1+L2+x(1,:)+x(2,:))+0.2;
%%
xplane=linspace(0,xmax,10);yplane=linspace(-0.1,0.1,10);
[Xplane,Yplane]=meshgrid(xplane,yplane);
[Yplane1,Zplane1]=meshgrid(yplane,yplane);
Zplane=-0.01*ones(size(Xplane));
Xplane1=zeros(size(Yplane));
Yplane2=0.1*ones(size(Yplane));
surf(Xplane,Yplane,Zplane,'FaceColor',[0 0.5 0],'EdgeColor','None');
hold on
surf(Xplane1,Yplane1,Zplane1,'FaceColor',[0 1 1],'EdgeColor','None');
surf(Xplane,Yplane2,Zplane1,'FaceColor',[0 1 1],'EdgeColor','None');
view(60, 45)
light
lighting phong
%%
plot3(xx,st,ct,'r')
plot3(xx2,st,ct,'r')
%%
surf(L1+x(1,i)+r*Xs,r*Ys,r*Zs,'FaceColor',[0 0 1],'EdgeColor','None');
surf(L1+L2+x(1,i)+x(2,i)+r*Xs,r*Ys,r*Zs,'FaceColor',[0 0 1],'EdgeColor','None');
axis([0 xmax -0.1 0.1 -0.01 0.1])
xlabel('x');ylabel('y');zlabel('z');
getframe;
hold off
end


% --- Executes on button press in pushbutton3.
function pushbutton3_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%extracting values of paremeter
m1=str2num(get(handles.edit1,'String'));
m2=str2num(get(handles.edit2,'String'));
L1=str2num(get(handles.edit3,'String'));
L2=str2num(get(handles.edit4,'String'));
k1=str2num(get(handles.edit5,'String'));
k2=str2num(get(handles.edit6,'String'));
x1=str2num(get(handles.edit13,'String'));
x2=str2num(get(handles.edit14,'String'));
v1=str2num(get(handles.edit15,'String'));
v2=str2num(get(handles.edit16,'String'));
T=str2num(get(handles.edit22,'String'));
h=str2num(get(handles.edit21,'String'));
g=str2num(get(handles.edit23,'String'));
mu1=str2num(get(handles.edit7,'String'));
mu2=str2num(get(handles.edit8,'String'));
alpha1=str2num(get(handles.edit9,'String'));
alpha2=str2num(get(handles.edit10,'String'));
beta1=str2num(get(handles.edit11,'String'));
beta2=str2num(get(handles.edit12,'String'));
F1=str2num(get(handles.edit17,'String'));
F2=str2num(get(handles.edit18,'String'));
w1=str2num(get(handles.edit19,'String'));
w2=str2num(get(handles.edit20,'String'));
I=str2num(get(handles.radiobutton1,'String'));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
y0=[x1;x2;v1;v2];%initial value
r=0.01;%radius for spring
t0=0;n=(T-t0)/h;t=linspace(t0,T,n);%assignining time values and calculating number os steps to be iterated
f=@(t,x) [x(3);x(4);-mu1*g-alpha1*x(3)+beta1*x(1)^3+beta2*(x(1)-x(2))^3-k1/m1*x(1)-k2/m1*(x(1)-x(2))+F1*cos(w1*t);-mu2*g-alpha2*x(4)+beta2*(x(2)-x(1))^3-k2/m2*(x(2)-x(1))+F2*cos(w2*t)];%defining the function
x=zeros(4,n);x(:,1)=y0;%defining the matrix for aproximations and initializing
[t1,Y]=ode45(f,t,y0);%calculating exact values
%defining RK-method
A=[0 0 0;0.5 0 0;-1 2 0];
c=[0 0.5 1];
b=[1/6 2/3 1/6];
xi=zeros(4,3);%RK-stages
xi(:,1)=y0;
for i=1:n-1
    if I==0
    xi(:,1)=y0;
    xi(:,2)=y0+h*A(2,1)*f(t(i)+c(1)*h,xi(:,1));
    xi(:,3)=y0+h*A(3,1)*f(t(i)+c(1)*h,xi(:,1))+h*A(3,2)*f(t(i)+c(2)*h,xi(:,2));
    y0=y0+h*(b(1)*f(t(i)+c(1)*h,xi(:,1))+b(2)*f(t(i)+c(2)*h,xi(:,2))+b(3)*f(t(i)+c(3)*h,xi(:,3)));
    x(:,i+1)=y0;
    else
        p=@(y) y-x(:,i)-h*f((t(i+1)+t(i))/2,(x(:,i)+y)/2);
       x(:,i+1)=fsolve(@(y)p(y),x(:,i));
    end
end
axes(handles.axes2)
 plot(t,x(1,:),'r-',t,Y(:,1),'b-');
 xlabel('Time(s)');
 ylabel('Elongation of spring 1');
 legend('Aproximation','Exact');
 title('Elongation of spring 1 vs time')
 axes(handles.axes3)
 plot(t,x(2,:),'r-',t,Y(:,2),'b-');
 xlabel('Time(s)');
 ylabel('Elongation of spring 2');
 legend('Aproximation','Exact');
 title('Elongation of spring 2 vs time');
for i=1:n-1
    axes(handles.axes1)
%%
tt = pi/2:pi/50:21*pi/2;
xx=linspace(0,L1+x(1,i),length(tt));
xx2=linspace(L1+x(1,i),L1+L2+x(1,i)+x(2,i),length(tt));
[Xs,Ys,Zs]=sphere(50);
st =0.01*sin(tt);
ct =0.01*cos(tt);
xmax=max(L1+L2+x(1,:)+x(2,:))+0.2;
%%
xplane=linspace(0,xmax,10);yplane=linspace(-0.1,0.1,10);
[Xplane,Yplane]=meshgrid(xplane,yplane);
[Yplane1,Zplane1]=meshgrid(yplane,yplane);
Zplane=-0.01*ones(size(Xplane));
Xplane1=zeros(size(Yplane));
Yplane2=0.1*ones(size(Yplane));
surf(Xplane,Yplane,Zplane,'FaceColor',[0 0.5 0],'EdgeColor','None');
hold on
surf(Xplane1,Yplane1,Zplane1,'FaceColor',[0 1 1],'EdgeColor','None');
surf(Xplane,Yplane2,Zplane1,'FaceColor',[0 1 1],'EdgeColor','None');
view(60, 45)
light
lighting phong
%%
plot3(xx,st,ct,'r')
plot3(xx2,st,ct,'r')
%%
surf(L1+x(1,i)+r*Xs,r*Ys,r*Zs,'FaceColor',[0 0 1],'EdgeColor','None');
surf(L1+L2+x(1,i)+x(2,i)+r*Xs,r*Ys,r*Zs,'FaceColor',[0 0 1],'EdgeColor','None');
axis([0 xmax -0.1 0.1 -0.01 0.1])
xlabel('x');ylabel('y');zlabel('z');
getframe;
hold off
end


% --- Executes on button press in radiobutton1.
function radiobutton1_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton1



function edit21_Callback(hObject, eventdata, handles)
% hObject    handle to edit21 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit21 as text
%        str2double(get(hObject,'String')) returns contents of edit21 as a double


% --- Executes during object creation, after setting all properties.
function edit21_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit21 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit22_Callback(hObject, eventdata, handles)
% hObject    handle to edit22 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit22 as text
%        str2double(get(hObject,'String')) returns contents of edit22 as a double


% --- Executes during object creation, after setting all properties.
function edit22_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit22 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit23_Callback(hObject, eventdata, handles)
% hObject    handle to edit23 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit23 as text
%        str2double(get(hObject,'String')) returns contents of edit23 as a double


% --- Executes during object creation, after setting all properties.
function edit23_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit23 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

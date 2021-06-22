function varargout = gui_spring3mass2(varargin)
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @gui_spring3mass2_OpeningFcn, ...
                   'gui_OutputFcn',  @gui_spring3mass2_OutputFcn, ...
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


function gui_spring3mass2_OpeningFcn(hObject, eventdata, handles, varargin)
global x1 x2 v1 v2 
global m1 m2 k1 k2 k3
global l1 l2 l3
global stop
stop=1;
num_masses=2;
%default values
m1=1; m2=1;
k1=1; k2=1; k3=1;
x1=0; x2=0;
v1=0; v2=0;
%length of springs
l1=1.1; l2=2.1; l3=1.1;
%set up initial figure of masses in equilibrium
axis([0 num_masses+l1+l2+l3 -1 1]);axis equal;axis off;set(gcf,'Color','w')
line([0 0],[-0.5 1],'Color','k'); hold on;
%line([num_masses+l1+l2+l3 num_masses+l1+l2+l3],[-0.5 1],'Color','k');
%line([0 num_masses+l1+l2+l3], [-0.5 -0.5],'Color','k');
%initial positions of masses
ypos=[-0.5 -0.5 0.5 0.5];
xpos=[l1+x1,1+l1+x1,1+l1+x1,l1+x1];
m1box=patch(xpos,ypos,'r'); 
xpos=[1+l1+l2+x2,2+l1+l2+x2,2+l1+l2+x2,1+l1+l2+x2];
m2box=patch(xpos,ypos,'r');
%initial position of springs
ne = 10; a = 1; ro = 0.1;
[xs,ys]=spring(0,0,l1+x1,0,ne,a,ro); spring1=plot(xs,ys,'LineWidth',2);
[xs,ys]=spring(1+l1+x1,0,1+l1+l2+x2,0,ne,a,ro); spring2=plot(xs,ys,'LineWidth',2);
[xs,ys]=spring(2+l1+l2+x2,0,2+l1+l2+l3,0,ne,a,ro);spring3=plot(xs,ys,'LineWidth',2);
%eigenvectors
A=[ [-(k1+k2), k2]/m1; [k2, -(k2+k3)]/m2 ];
[eigenvectors, eigenvalues]=eig(A);
[~,ix]=sort(diag(eigenvalues));
eigenvectors=eigenvectors(:,ix); 
%data to share
handles.m1box=m1box;handles.m2box=m2box;
handles.spring1=spring1;handles.spring2=spring2;handles.spring3=spring3;
handles.l1=l1;handles.l2=l2;handles.l3=l3;
handles.eigenvectors=eigenvectors;


handles.output = hObject;

% Update handles
guidata(hObject, handles);


function varargout = gui_spring3mass2_OutputFcn(hObject,eventdata,handles) 

varargout{1} = handles.output;

%runs spring
function run(handles)
global x1 x2 v1 v2
global m1 m2 k1 k2 k3
global stop

N=64; T=(2*pi); dt=T/N; dt_pause=0.016;

%shared data 
m1box=handles.m1box;m2box=handles.m2box;
spring1=handles.spring1;spring2=handles.spring2;spring3=handles.spring3;
l1=handles.l1;l2=handles.l2;l3=handles.l3;

stop=0;


%curve = animatedline('Color','r','Marker','o');
%set(gca,'XLim',[0 2*pi], 'YLim',[-2 2]);
%grid on;

lista = [];
lista2 = [];
x = 0:0.05:2*pi;

for i=1:intmax
    t_loopstart=tic();
    if stop
        break;
    end
    %fprintf('%g %f %f %f \n',i,x1,x2,x3);
    [t,y] = ode45(@(t,y) masses(t,y,m1,m2,k1,k2,k3),...
                                               [0,dt],[x1,v1,x2,v2]);
 
    x1=y(end,1); x2=y(end,3); 
    v1=y(end,2); v2=y(end,4);
    xpos=[l1+x1,1+l1+x1,1+l1+x1,l1+x1];
    set(m1box,'xdata',xpos);
    xpos=[1+l1+l2+x2,2+l1+l2+x2,2+l1+l2+x2,1+l1+l2+x2];
    set(m2box,'xdata',xpos); 
    [xs,ys]=spring(0,0,l1+x1,0); set(spring1,'xdata',xs,'ydata',ys);
    [xs,ys]=spring(1+l1+x1,0,1+l1+l2+x2,0); set(spring2,'xdata',xs,'ydata',ys);
    [xs,ys]=spring(2+l1+l2+x2,0,2+l1+l2+l3,0); set(spring3,'xdata',xs,'ydata',ys);
    el_time=toc(t_loopstart);
    %disp(['Elapse time : ',num2str(el_time),' second']);
    
    
    x(i+1) = x(i)+0.05;
    lista(end+1) = x1;
    lista2(end+1) = x2;
    
    if i<=200
        %plot(x(i),x1,'g','Parent',handles.graph);
        plot(x(1:i),lista(1:i),'g','Parent',handles.graph)
        hold(handles.graph, 'on')
        %plot(x(i),x2,'r','Parent',handles.graph);
        plot(x(1:i),lista2(1:i),'r','Parent',handles.graph);
        hold(handles.graph, 'on')
        
    else
        hold(handles.graph, 'off')
        %plot(x(i),x1,'g','Parent',handles.graph); 
        plot(x(i-200:i),lista(i-200:i),'g','Parent',handles.graph);
        hold(handles.graph, 'on')
        %plot(x(i),x2,'r','Parent',handles.graph); 
        plot(x(i-200:i),lista2(i-200:i),'r','Parent',handles.graph);
        hold(handles.graph, 'on')
    end
    
    %if i<=200
    %    plot(x(i),x2,'r','Parent', handles.graph);
    %    hold on;
    %    plot(x(1:i),lista2(1:i),'r','Parent', handles.graph);
    %else
    %    plot(x(i),x2,'r','Parent', handles.graph); 
    %    hold on;
    %    plot(x(i-200:i),lista2(i-200:i),'r','Parent', handles.graph);
    %end
    
    pause(dt_pause-el_time);
end

function dy=masses(~,y,m1,m2,k1,k2,k3)
dy=zeros(4,1);
dy(1)=y(2);
dy(2)=(-(k1+k2)*y(1) + k2*y(3))/m1;
dy(3)=y(4);
dy(4)=(k2*y(1) - (k2+k3)*y(3))/m2;

function m1_Callback(hObject, eventdata, handles)

global m1 m2 k1 k2 k3 
global stop
stop=1;
m1temp=str2double(get(hObject, 'String'));
if m1temp <= 0
    msgbox('Parameter out-of-range.  Choose m1 > 0.');
    set(hObject, 'String', num2str(m1));
    return;
else
    m1=m1temp;
end

A=[ [-(k1+k2), k2]/m1; [k2, -(k2+k3)]/m2 ];
[eigenvectors, eigenvalues]=eig(A);
[~,ix]=sort(diag(eigenvalues));
eigenvectors=eigenvectors(:,ix); %sort from higher frequency to lower
handles.eigenvectors=eigenvectors;
guidata(hObject,handles)

function m1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to m1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function m2_Callback(hObject, eventdata, handles)
% hObject    handle to m2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of m2 as text
%      str2double(get(hObject,'String')) returns contents of m2 as a double
global m1 m2 k1 k2 k3
global stop
stop=1;

m2temp=str2double(get(hObject, 'String'));
if m2temp <= 0
    msgbox('Parameter out-of-range.  Choose m2 > 0.');
    set(hObject, 'String', num2str(m2));
    return;
else
    m2=m2temp;
end

A=[ [-(k1+k2), k2]/m1; [k2, -(k2+k3)]/m2 ];
[eigenvectors, eigenvalues]=eig(A);
[~,ix]=sort(diag(eigenvalues));
eigenvectors=eigenvectors(:,ix); %sort from higher frequency to lower
handles.eigenvectors=eigenvectors;
guidata(hObject,handles)


function m2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to m2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function k1_Callback(hObject, eventdata, handles)
% hObject    handle to k1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of k1 as text
%      str2double(get(hObject,'String')) returns contents of k1 as a double
global m1 m2 k1 k2 k3
global stop
stop=1;

k1temp=str2double(get(hObject, 'String'));
if k1temp <= 0
    msgbox('Parameter out-of-range.  Choose k1 > 0.');
    set(hObject, 'String', num2str(k1));
    return;
else
    k1=k1temp;
end

A=[ [-(k1+k2), k2]/m1; [k2, -(k2+k3)]/m2 ];
[eigenvectors, eigenvalues]=eig(A);
[~,ix]=sort(diag(eigenvalues));
eigenvectors=eigenvectors(:,ix); %sort from higher frequency to lower
handles.eigenvectors=eigenvectors;
guidata(hObject,handles)


function k1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to k1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function k2_Callback(hObject, eventdata, handles)
% hObject    handle to k2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of k2 as text
%      str2double(get(hObject,'String')) returns contents of k2 as a double
global m1 m2 k1 k2 k3
global stop
stop=1;

k2temp=str2double(get(hObject, 'String'));
if k2temp <= 0
    msgbox('Parameter out-of-range.  Choose k2 > 0.');
    set(hObject, 'String', num2str(k2));
    return;
else
    k2=k2temp;
end
A=[ [-(k1+k2), k2]/m1; [k2, -(k2+k3)]/m2 ];
[eigenvectors, eigenvalues]=eig(A);
[~,ix]=sort(diag(eigenvalues));
eigenvectors=eigenvectors(:,ix); %sort from higher frequency to lower
handles.eigenvectors=eigenvectors;
guidata(hObject,handles)


function k2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to k2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function k3_Callback(hObject, eventdata, handles)
% hObject    handle to k3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of k3 as text
%      str2double(get(hObject,'String')) returns contents of k3 as a double
global m1 m2 k1 k2 k3
global stop
stop=1;

k3temp=str2double(get(hObject, 'String'));
if k3temp <= 0
    msgbox('Parameter out-of-range.  Choose k3 > 0.');
    set(hObject, 'String', num2str(k3));
    return;
else
    k3=k3temp;
end

A=[ [-(k1+k2), k2]/m1; [k2, -(k2+k3)]/m2 ];
[eigenvectors, eigenvalues]=eig(A);
[~,ix]=sort(diag(eigenvalues));
eigenvectors=eigenvectors(:,ix); %sort from higher frequency to lower
handles.eigenvectors=eigenvectors;
guidata(hObject,handles)


function k3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to k3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function first_eigenmode_Callback(hObject, eventdata, handles)
% hObject    handle to first_eigenmode (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global x1 x2 v1 v2
global stop

x1=handles.eigenvectors(1,1);
x2=handles.eigenvectors(2,1);
%v1=0; v2=0;

if stop; run(handles); end;


function second_eigenmode_Callback(hObject, eventdata, handles)
% hObject    handle to second_eigenmode (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global x1 x2 v1 v2
global stop

x1=handles.eigenvectors(1,2);
x2=handles.eigenvectors(2,2);

%v1=0; v2=0;

if stop; run(handles); end;



function random_Callback(hObject, eventdata, handles)
% hObject    handle to random (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%if running, stop it
global x1 x2  v1 v2
global stop

eigenvectors=handles.eigenvectors;
coefficients=2*(rand(2,1)-0.5);
xvec=eigenvectors*coefficients;
x1=xvec(1); x2=xvec(2);
v1=0; v2=0;

if stop; run(handles); end;


function quit_Callback(hObject, eventdata, handles)
% hObject    handle to quit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global stop
stop=1;
close all;
closereq();



function edit7_Callback(hObject, eventdata, handles)
% hObject    handle to edit7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit7 as text
%        str2double(get(hObject,'String')) returns contents of edit7 as a double
global l1 l2 l3
global x1 x2


l1temp=str2double(get(hObject, 'String'));
if l1temp <= 0
    msgbox('Parameter out-of-range.  Choose l1 > 0.');
    set(hObject, 'String', num2str(l1));
    return;
else
    l1=l1temp;

    
first_pos = [0,6.3,6.3,0];
second_pos = [-0.5,-0.5,1,1];
white_box = patch(first_pos, second_pos,'w', EdgeColor = 'w');
%initial positions of masses
ypos=[-0.5 -0.5 0.5 0.5];
xpos=[l1+x1,1+l1+x1,1+l1+x1,l1+x1];
m1box=patch(xpos,ypos,'r'); 
xpos=[1+l1+l2+x2,2+l1+l2+x2,2+l1+l2+x2,1+l1+l2+x2];
m2box=patch(xpos,ypos,'r');
%initial position of springs
ne = 10; a = 1; ro = 0.1;
[xs,ys]=spring(0,0,l1+x1,0,ne,a,ro); spring1=plot(xs,ys,'LineWidth',2);
[xs,ys]=spring(1+l1+x1,0,1+l1+l2+x2,0,ne,a,ro); spring2=plot(xs,ys,'LineWidth',2);
[xs,ys]=spring(2+l1+l2+x2,0,2+l1+l2+l3,0,ne,a,ro);spring3=plot(xs,ys,'LineWidth',2);

handles.m1box=m1box;handles.m2box=m2box;
handles.spring1=spring1;handles.spring2=spring2;handles.spring3=spring3;
handles.l1=l1;handles.l2=l2;handles.l3=l3;

guidata(hObject, handles);
end


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
global l1 l2 l3
global x1 x2

l2temp=str2double(get(hObject, 'String'));
if l2temp <= 0
    msgbox('Parameter out-of-range.  Choose l2 > 0.');
    set(hObject, 'String', num2str(l2));
    return;
else
    l2=l2temp;

first_pos = [0,6.3,6.3,0];
second_pos = [-0.5,-0.5,1,1];
white_box = patch(first_pos, second_pos,'w', EdgeColor = 'w');
%initial positions of masses
ypos=[-0.5 -0.5 0.5 0.5];
xpos=[l1+x1,1+l1+x1,1+l1+x1,l1+x1];
m1box=patch(xpos,ypos,'r'); 
xpos=[1+l1+l2+x2,2+l1+l2+x2,2+l1+l2+x2,1+l1+l2+x2];
m2box=patch(xpos,ypos,'r');
%initial position of springs
ne = 10; a = 1; ro = 0.1;
[xs,ys]=spring(0,0,l1+x1,0,ne,a,ro); spring1=plot(xs,ys,'LineWidth',2);
[xs,ys]=spring(1+l1+x1,0,1+l1+l2+x2,0,ne,a,ro); spring2=plot(xs,ys,'LineWidth',2);
[xs,ys]=spring(2+l1+l2+x2,0,2+l1+l2+l3,0,ne,a,ro);spring3=plot(xs,ys,'LineWidth',2);

handles.m1box=m1box;handles.m2box=m2box;
handles.spring1=spring1;handles.spring2=spring2;handles.spring3=spring3;
handles.l1=l1;handles.l2=l2;handles.l3=l3;

guidata(hObject, handles);
end


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

global l1 l2 l3
global x1 x2


l3temp=str2double(get(hObject, 'String'));
if l3temp <= 0
    msgbox('Parameter out-of-range.  Choose l3 > 0.');
    set(hObject, 'String', num2str(l3));
    return;
else
    l3=l3temp;

first_pos = [0,6.3,6.3,0];
second_pos = [-0.5,-0.5,1,1];
white_box = patch(first_pos, second_pos,'w', EdgeColor = 'w');
%initial positions of masses
ypos=[-0.5 -0.5 0.5 0.5];
xpos=[l1+x1,1+l1+x1,1+l1+x1,l1+x1];
m1box=patch(xpos,ypos,'r'); 
xpos=[1+l1+l2+x2,2+l1+l2+x2,2+l1+l2+x2,1+l1+l2+x2];
m2box=patch(xpos,ypos,'r');
%initial position of springs
ne = 10; a = 1; ro = 0.1;
[xs,ys]=spring(0,0,l1+x1,0,ne,a,ro); spring1=plot(xs,ys,'LineWidth',2);
[xs,ys]=spring(1+l1+x1,0,1+l1+l2+x2,0,ne,a,ro); spring2=plot(xs,ys,'LineWidth',2);
[xs,ys]=spring(2+l1+l2+x2,0,2+l1+l2+l3,0,ne,a,ro);spring3=plot(xs,ys,'LineWidth',2);

handles.m1box=m1box;handles.m2box=m2box;
handles.spring1=spring1;handles.spring2=spring2;handles.spring3=spring3;
handles.l1=l1;handles.l2=l2;handles.l3=l3;

guidata(hObject, handles);
end



function edit9_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit9 (see GCBO)
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
global v1

v1temp=str2double(get(hObject, 'String'));
if v1temp <= 0
    msgbox('Parameter out-of-range.  Choose l3 > 0.');
    set(hObject, 'String', num2str(l3));
    return;
else
    v1=v1temp;
end


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
global v2

v2temp=str2double(get(hObject, 'String'));
if v2temp <= 0
    msgbox('Parameter out-of-range.  Choose l3 > 0.');
    set(hObject, 'String', num2str(l3));
    return;
else
    v2=v2temp;
end


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

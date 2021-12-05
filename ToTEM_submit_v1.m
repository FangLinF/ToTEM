function varargout = ToTEM_submit_v1(varargin)
% TOTEM_SUBMIT_V1 MATLAB code for ToTEM_submit_v1.fig

% Last Modified by GUIDE v2.5 27-May-2021 11:47:40

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @ToTEM_submit_v1_OpeningFcn, ...
                   'gui_OutputFcn',  @ToTEM_submit_v1_OutputFcn, ...
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


% --- Executes just before ToTEM_submit_v1 is made visible.
function ToTEM_submit_v1_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to ToTEM_submit_v1 (see VARARGIN)

% Choose default command line output for ToTEM_submit_v1
handles.output = hObject;

handles.exepath=cd; %��סexe�ļ����ڵ��ļ���
handles.execd=cd;  %��¼exe�ļ����ڵ�Ŀ¼
handles.tempcd=cd;  %������¼��һ���򿪵�Ŀ¼λ�ã������ļ���д
handles.saveresult = cd; % remember the current dir to save the result

set(handles.text21, 'visible', 'off');   % aper��mrad��1/nm��λ
set(handles.text125, 'visible', 'on');
set(handles.text24, 'visible', 'on');   % convergence & aperture��mrad��1/nm��λ
set(handles.text134, 'visible', 'off');
    
axis(handles.axes1); axis off;axis(handles.axes2); axis off;


% %x��y������Ҫ������չ4��
handles.outside_ext=4;   
%����Ǽ���ÿ��ԭ�ӵ��Ƴ�������ÿ��ԭ�ӵļ���
handles.TOTALELEment=10;
%�����񶯵Ĵ���
handles.vibration=30;
%��ʾԭ��ʱ��ԭ�ӵĳߴ�
handles.dis_atomsize = 1;

%һ���ܹ�������ٸ�probe�Ĵ�����
handles.mynode=70;

%STEM����
set(handles.radiobutton11,'value',1)
% Update handles structure

handles.conv_source=0; 
handles.conv_sampling=1; %�Ѿ���ĵ�Դ�ߴ磬�Լ��Ŵ��ʳ�����0��1�ĳ�ֵ�������κβ���

set(handles.GPURB,'value',1);  %Ĭ��GPU����


%��������ͼ��һЩ��Ϣ���������趨��ֵ
%�ļ�·��
handles.batchs.PathName = cd;
handles.batchs.wavesx=320;  handles.batchs.wavesy=handles.batchs.wavesx;
handles.batchs.showorrun = 'R';
handles.batchs.totalnumber =500;
handles.batchs.processnum = 2;
handles.batchs.GPUBatch = 5;
handles.batchs.top_left = [101,51];
handles.batchs.width_heigh = [50,70];

guidata(hObject, handles);



% --- Outputs from this function are returned to the command line.
function varargout = ToTEM_submit_v1_OutputFcn(hObject, eventdata, handles) 
varargout{1} = handles.output;



% --------------------------------------------------------------------
function Untitled_7_Callback(hObject, eventdata, handles)
function Untitled_4_Callback(hObject, eventdata, handles)
function Untitled_10_Callback(hObject, eventdata, handles)
function Untitled_1_Callback(hObject, eventdata, handles)

% --------------------------------------------------------------------
function Untitled_5_Callback(hObject, eventdata, handles)
cd(handles.tempcd)  %������һ�δ򿪵�Ŀ¼
[FileName,PathName]=uigetfile({'*.pdb','Program Database File (*.pdb)'});%���ļ�
if length(FileName)==1 & FileName(1)==0;  %���ѡ��cancalû���ļ��Ļ���������ֹ
    cd(handles.execd)
    return;
end
set(handles.edit1,'string',FileName);

cd(handles.execd);  %�ص��ɵ��ļ���
handles.tempcd=PathName;  %��¼��һ���򿪵��ļ���

%��ȡpdb����cif�ľ�������
if sum(FileName(end-2:end)=='pdb')==3  %��ʾ��pdb��ʽ
    xxx=readpdb(strcat(PathName, FileName));
else
    disp('Can not read this file');
    return;  
end
handles.x=[];
handles.x(:,1:4)=xxx(:,1:4);  %element��coordinates
handles.x(:,5) = 0;  %��̬
handles.x(:,6) = 0;  %������
handles.x(:,7) = xxx(:,6); %dw %will input from interface
handles.x(:,8) = xxx(:,5); %occupy

%input DW factors
    z=unique(handles.x(:,1));
    for i=1:length(z)
        prompt={sprintf('DW factor for the %.0f Atom ',z(i))};
        dlg_title = 'Input';
        num_lines = 1;
        defaultans = {'0'};
        DW=inputdlg(prompt,dlg_title,num_lines,defaultans);
        handles.x(find(handles.x(:,1)==z(i)),7)=str2num(cell2mat(DW));
    end
    temp = handles.x;
%     xlswrite(strcat(PathName, FileName(1:end-4)),temp);
%     disp(strcat('Parameters of atoms are saved in the file', strcat(PathName, FileName(1:end-4)),'.xlsx'));

%handles.x is a variable including all information of atoms, the atoms'
%coordinators is in the range of [0, MaxValue]. Specially, handles.x will
%not be changed in following codes.------begin 202012290930
atompos=handles.x(:,2:4);  % ��λ�䰣~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
atompos(:,3) = atompos(:,3)-min(atompos(:,3)); %+1.5*eachthick;  % ԭ�Ӹ߶ȣ���Сֵ��Ϊ0��������һ����ղ��ȣ�����ֲ�
                                                            % ���ұ�֤����λ�ڵ�0�㡣��Ϊ����������λ�ڵڶ��㣬������ȡ��floor��������
                                                            % ��������3,0-3-6-9������1.5*eachslice��ֵΪ4.5�����ڵ�һ�����м�

atompos(:,1) = atompos(:,1)-min(atompos(:,1)); %+handles.extendis;
atompos(:,2) = atompos(:,2)-min(atompos(:,2)); %+handles.extendis;

handles.atompos = atompos;  %�����޸�һ��ԭ�ӵ����꣬������Ϣ��Ҫ�鿴handles.x�����ȡ�Ľ����
                %ֻ��handles.x�е�������Ҫ�ı�һ�£�ֻ�������ò��ʱ��
handles.x(:,2:4)=atompos;  %���¸�ԭ�ӵ����꣬��֤���Ǵ���0�ģ����б߽�
%%%%%%%%%%%%%%%%%%%%%%%-----------------end 202012290930

% the atoms recoded in handles.x_new is just for draw. 
% handles.x_new will be reedit according to the view point and crystal
% rotation. 
%Also see the function of 'Rotation & view' and 'Projected along view'
%---begin 202012290935
handles.x_new=handles.x;  %����תǰ�Ľṹ��¼������handles.x_new��¼������ת��
Drawsupercell(hObject, eventdata, handles)
%%%%%%%%%%%%%%%%%%%%%%-----------------end 202012290935
guidata(hObject, handles);


function Drawsupercell(hObject, eventdata, handles)
%����ԭ�ӽṹ draw atoms of crystal in axis1 and axis2 from handles.x_new
%��Ҫ�����Ƚ϶����Ϣ�������������pdb��Ϣ����,��һ����Ԫ���������������������꣬������ռ���ʣ������м�װ�Ǽ�̬
atompos=handles.x_new(:,2:4);
axis(handles.axes1); hold off
cla(handles.axes1);   %���ԭ��ͼ��
cla(handles.axes2);   %���ԭ��ͼ��

%��ͼ��ʾһ�½ṹ
axes(handles.axes1);   
hold on;
view([0 0 -1])

maxx=max(atompos(:,1));
maxy=max(atompos(:,2));
maxz=max(atompos(:,3));
mm=max([maxx, maxy, maxz]);

hold on; line([0 0], [0 0], [0, mm],'color','b')
line([0 0], [0 mm], [0, 0],'color','g')
line([0 mm], [0 0], [0, 0],'color','r')
line([0 0], [0 0], [0, -mm/10],'color','b','linestyle','--')
line([0 0], [0 -mm/10], [0, 0],'color','g','linestyle','--')
line([0 -mm/10], [0 0], [0, 0],'color','r','linestyle','--')

item_ele=unique(handles.x(:,1));  %�����м���ԭ��
for i=1:length(item_ele)  %Ԫ�ص����� there are many elements
    r=1-mod(item_ele(i),5)*0.25;    %����ԭ�����������ø�ԭ�����е���ɫ��r,g,bΪ����ɫ��ֵ
    g=1-mod(floor(item_ele(i)/5),5)*0.25;
    b=1-mod(floor(item_ele(i)/25),5)*0.25;
    %���л�ͼ
    axes(handles.axes1);  %��ԭ�� draw atoms
    plot3(atompos(find(handles.x(:,1)==item_ele(i)),1),...
        atompos(handles.x(:,1)==item_ele(i),2), ...
        atompos(find(handles.x(:,1)==item_ele(i)),3),'o',...  %��㣬��״Ϊ��o��
         'MarkerEdgeColor',[0 0 0],...    %���ñ�Ե��ɫΪ��ɫ
         'MarkerSize',handles.dis_atomsize.*(20-log2(2)),...  %����ԭ�ӵĴ�С,������ʾcell��Ŀ�����Ӵ�С���С
         'MarkerFaceColor',[r g b]);    %����ԭ�ӵ���ɫ
     axes(handles.axes2);  hold on;%��ͼ��ʾһ��ͼ�� draw atoms' label
     plot(0.1*i,0.7,'o',...     %��ʾԭ��
     'MarkerEdgeColor',[0 0 0],...
     'MarkerSize',15,...
     'MarkerFaceColor',[r g b]);
     text(0.1*i-0.01,0.3,num2str(item_ele(i)));  %��ԭ���·���ʾ��Ӧ��ԭ��������
end
axes(handles.axes1);  %��ԭ��
axis equal
xlabel('x'); ylabel('y'); zlabel('z');
axes(handles.axes2);
box off;axis off;   %����ʾ�߿��������
xlim([0 1]);ylim([0 1]);  %����x����y���ȡֵ��Χ
guidata(hObject, handles);
   


% --------------------------------------------------------------------
function Untitled_6_Callback(hObject, eventdata, handles)
cd(handles.tempcd)  %������һ�δ򿪵�Ŀ¼
[FileName,PathName]=uigetfile({'*.xls','Excel File (*.xls)';'*.xlsx','Excel File (*.xlsx)'});%���ļ�
if length(FileName)==1 & FileName(1)==0;  %���ѡ��cancalû���ļ��Ļ���������ֹ
    cd(handles.execd)
    return;
end
set(handles.edit1,'string',FileName);

cd(handles.execd);  %�ص��ɵ��ļ���
handles.tempcd=PathName;  %��¼��һ���򿪵��ļ���

%��ȡpdb����cif�ľ�������
if sum(FileName(end-2:end)=='xls')==3 || sum(FileName(end-3:end)=='xlsx')==4 %��ʾ��pdb��ʽ
    handles.x=xlsread(strcat(PathName, FileName));
else
    disp('Cannot read this file!');
    return;
end

atompos=handles.x(:,2:4);  % ��λ�䰣~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
atompos(:,3) = atompos(:,3)-min(atompos(:,3)); %+1.5*eachthick;  % ԭ�Ӹ߶ȣ���Сֵ��Ϊ0��������һ����ղ��ȣ�����ֲ�
                                                            % ���ұ�֤����λ�ڵ�0�㡣��Ϊ����������λ�ڵڶ��㣬������ȡ��floor��������
                                                            % ��������3,0-3-6-9������1.5*eachslice��ֵΪ4.5�����ڵ�һ�����м�

atompos(:,1) = atompos(:,1)-min(atompos(:,1)); %+handles.extendis;
atompos(:,2) = atompos(:,2)-min(atompos(:,2)); %+handles.extendis;

handles.atompos = atompos;  %�����޸�һ��ԭ�ӵ����꣬������Ϣ��Ҫ�鿴handles.x�����ȡ�Ľ����
                %ֻ��handles.x�е�������Ҫ�ı�һ�£�ֻ�������ò��ʱ��
handles.x(:,2:4)=atompos;  %���¸�ԭ�ӵ����꣬��֤���Ǵ���0�ģ����б߽�
%%%%%%%%%%%%%%%%%%%%%%%-----------------end 202012290930

% the atoms recoded in handles.x_new is just for draw. 
% handles.x_new will be reedit according to the view point and crystal
% rotation. 
%Also see the function of 'Rotation & view' and 'Projected along view'
%---begin 202012290935
handles.x_new=handles.x;  %����תǰ�Ľṹ��¼������handles.x_new��¼������ת��
Drawsupercell(hObject, eventdata, handles)
guidata(hObject, handles);


function radiobutton2_Callback(hObject, eventdata, handles)
% if get(handles.radiobutton2,'value')==1
%     set(handles.radiobutton4, 'value', 0);
% end
guidata(hObject, handles);




function radiobutton4_Callback(hObject, eventdata, handles)
% if get(handles.radiobutton4,'value')==1
%     set(handles.radiobutton2, 'value', 0);
% end
guidata(hObject, handles);



% --- Executes on button press in checkbox2.
function checkbox1_Callback(hObject, eventdata, handles)
if get(handles.popupmenu5,'value') == 3
    set(handles.popupmenu5,'value',2);  %�����Ҫ�ֲ���㣬��������ת��lobato��ϵ������Ҫǿ�л�Ϊpeng��
end
guidata(hObject, handles);


function edit4_Callback(hObject, eventdata, handles)
% --- Executes during object creation, after setting all properties.
function edit4_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit20_Callback(hObject, eventdata, handles)
% --- Executes during object creation, after setting all properties.
function edit20_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit19_Callback(hObject, eventdata, handles)
% --- Executes during object creation, after setting all properties.
function edit19_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit18_Callback(hObject, eventdata, handles)
function edit18_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit17_Callback(hObject, eventdata, handles)
function edit17_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_Vol_Callback(hObject, eventdata, handles)
function edit_Vol_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_Con_Callback(hObject, eventdata, handles)
function edit_Con_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_Ape_Callback(hObject, eventdata, handles)
function edit_Ape_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_Spr_Callback(hObject, eventdata, handles)
function edit_Spr_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_Tilt_1_Callback(hObject, eventdata, handles)
function edit_Tilt_1_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_Tilt_2_Callback(hObject, eventdata, handles)
function edit_Tilt_2_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_Focus_Callback(hObject, eventdata, handles)
function edit_Focus_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_A1_1_Callback(hObject, eventdata, handles)
function edit_A1_1_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_A1_2_Callback(hObject, eventdata, handles)
function edit_A1_2_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_A2_1_Callback(hObject, eventdata, handles)
function edit_A2_1_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_A2_2_Callback(hObject, eventdata, handles)
function edit_A2_2_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_B2_1_Callback(hObject, eventdata, handles)
function edit_B2_1_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_B2_2_Callback(hObject, eventdata, handles)
function edit_B2_2_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_Cs_Callback(hObject, eventdata, handles)
h_ba=(6.6255916e-34)/2/pi;
h=6.6255916e-34;
e=1.602102e-19;
me=9.1093807e-31;
c=2.9979251e+8;
%����scherzer����
vol=str2num(get(handles.edit_Vol,'string'));
lambda=1.0e+9*h.*c./sqrt(e*vol*1000*(2*me*c*c+e*vol*1000));  %���㲨������λnm
str=num2str(-1.2*sqrt( 1000*str2num(get(handles.edit_Cs, 'string'))*lambda));
disp(strcat('Scherzer focus is:',str, 'nm'));

% --- Executes during object creation, after setting all properties.
function edit_Cs_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_A3_1_Callback(hObject, eventdata, handles)
function edit_A3_1_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_A3_2_Callback(hObject, eventdata, handles)
function edit_A3_2_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_S3_1_Callback(hObject, eventdata, handles)
function edit_S3_1_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_S3_2_Callback(hObject, eventdata, handles)
function edit_S3_2_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_A4_1_Callback(hObject, eventdata, handles)
function edit_A4_1_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_A4_2_Callback(hObject, eventdata, handles)
function edit_A4_2_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_D4_1_Callback(hObject, eventdata, handles)
function edit_D4_1_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_D4_2_Callback(hObject, eventdata, handles)
function edit_D4_2_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_B4_1_Callback(hObject, eventdata, handles)
function edit_B4_1_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_B4_2_Callback(hObject, eventdata, handles)
function edit_B4_2_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_C5_Callback(hObject, eventdata, handles)
function edit_C5_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_A5_1_Callback(hObject, eventdata, handles)
function edit_A5_1_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_A5_2_Callback(hObject, eventdata, handles)
function edit_A5_2_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_D5_1_Callback(hObject, eventdata, handles)
function edit_D5_1_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_D5_2_Callback(hObject, eventdata, handles)
function edit_D5_2_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_S5_1_Callback(hObject, eventdata, handles)
function edit_S5_1_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_S5_2_Callback(hObject, eventdata, handles)
function edit_S5_2_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit75_Callback(hObject, eventdata, handles)
function edit75_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit77_Callback(hObject, eventdata, handles)
function edit77_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%----------------------------------------
function resul=Gaussian_focal(mtotal,rms,lambda,gmax);   %�����˹�ֲ���ֵ
M=(mtotal-1)/2;
if M==0
    resul.delta_yita=0;  %�뽹ƫ���Ƕ���
    resul.gfsf=1;      %��˹ȡֵ�Ƕ���
else
    resul.delta_yita=energyspread(lambda, mtotal, rms, gmax);  %��Ϊ����ͼ�ɶ���ͼ���ӣ�ȷ���뽹���ľ�����ֵ
    resul.gfsf=gaussian_focal(mtotal,rms,resul.delta_yita);   %ȷ�������뽹�����£���˹��ɢ�İٷֱȱ��� 
end


function f_delta=gaussian_focal(mtotal,delta,delta_yita);%the calculation of a Gaussian focal spread function
M=(mtotal-1)/2;
if M==0
    f_delta=1;
else
   f_delta=zeros(1,2*M+1);
   a=zeros(1,2*M+1);
   a=(-M:M)*delta_yita;
   f_delta=delta_yita*exp(-a.*a/(2*delta*delta))/(sqrt(2*pi)*delta);
end


function delta_yita=energyspread(lambda, mtotal, rms, gmax)   %�������ŵ�deltayita��ֵ�����뽹ƫ���Ƕ���
M=(mtotal-1)/2;
%yita ----- yita is the focal_step for rms of gaussian function
%pre  ----- the precision for the yita
pre=0.1; %unit nm
fieldmin=0.5;fieldmax=2.0;  %----- search yita from rms-field,unit is nm
uu=0:0.01:gmax;
fixvalue=exp(-0.5*(pi*rms*lambda).^2*uu.^4);
tempdelta=fieldmin.*rms:pre:fieldmax.*rms;
for i=1:length(tempdelta)
    gfsf=gaussian_focal(mtotal,rms,tempdelta(i));
    f_delta_value(i)=max( abs( fixvalue-real( ( gfsf*exp( -sqrt( -1 ) *pi*lambda*( ( -M:M )' )*tempdelta(i)*uu.^2) ) ) ) );   %���֮��Ĳ�ֵ,ע������������������¼����󣬳�������ȷ��
end
delta_yita=tempdelta(find(f_delta_value==min(f_delta_value)));
pp=1;



function tccpara=readtccfromhandles_newSTEM(hObject, handles, lambda)  %�����Ǹ���TEM�Ĳ�������ȡ��
tccpara.lambda=lambda;
%if get(handles.checkbox_polar,'Value')==1  %�ض��Ǽ�����ϵ
   %��ȡ˵����%�ο�Phil. Trans. R. Soc. A 2009, vol
%367:3755-3771;����Ȧ������ϵ����UM  1998 72 PP109-119�Լ�UM 1996 64 249-264�У�B2 S3 D4 B4
%S5 R5(D5)����û��ϵ������ˣ���ȡ���ݺ�Ӧ��Ҫ�ֱ���� 3  4  5  5  6  6

tccpara.focus=str2num(get(handles.edit_Focus,'String'));
   
   %A1����ֵ���䣻A1�ĽǶ��ǵ羵��ʾ����ֵ��1/2 
tccpara.A1=str2num(get(handles.edit_A1_1,'String')); tccpara.phiA1=-str2num(get(handles.edit_A1_2,'String'));
   
   %A2����ֵ���䣻A2�ĽǶ��ǵ羵��ʾ����ֵ��1/3
   tccpara.A2=str2num(get(handles.edit_A2_1,'String')); tccpara.phiA2=-str2num(get(handles.edit_A2_2,'String'));
  
   %B2����ֵ���䣻B2�ĽǶ���羵��ʾ��ȡ����
   tccpara.B2=str2num(get(handles.edit_B2_1,'String')); %�����������ʱ���rew�����3��������rew����Ϊ2000�����������6000����Ϊ�����е�ϵ�������1/3
   tccpara.phiB2=-str2num(get(handles.edit_B2_2,'String'));  
  
   %Cs�ĵ�λ��um��nm
   tccpara.Cs=str2num(get(handles.edit_Cs,'String'))*(10^3); %UNIT from um to nm
   
   %A3�ĵ�λ��um��nm��A3�ĽǶ���羵��ʾ�ĳ���4
   tccpara.A3=str2num(get(handles.edit_A3_1,'String'))*(10^3); %UNIT from um to nm
   tccpara.phiA3=-str2num(get(handles.edit_A3_2,'String'));
   
   %S3�ĵ�λ��um��nm��S3�ĽǶ���羵��ʾ��ȡ��������2
   tccpara.S3=str2num(get(handles.edit_S3_1,'String'))*(10^3); 
   tccpara.phiS3=-str2num(get(handles.edit_S3_2,'String'));
   
   %A4�ĵ�λ��um��nm��A4�ĽǶ���羵��ʾ�ĳ���5
   tccpara.A4=str2num(get(handles.edit_A4_1,'String'))*(10^3);
   tccpara.phiA4=-str2num(get(handles.edit_A4_2,'String'));
   
   %D4�ĵ�λ��um��nm��D4�ĽǶ���羵��ʾ��ȡ������3
   tccpara.D4=str2num(get(handles.edit_D4_1,'String'))*(10^3); 
   tccpara.phiD4=-str2num(get(handles.edit_D4_2,'String'));
   
   %B4�ĵ�λ��um��nm��B4�ĽǶ���羵��ʾ��һ��   %��Ҫ���о��ĽǶȹ�ϵ������
   tccpara.B4=str2num(get(handles.edit_B4_1,'String'))*(10^3); 
   tccpara.phiB4=-str2num(get(handles.edit_B4_2,'String'));
   
    %A5�ĵ�λ��mm��nm��A4�ĽǶ���羵��ʾ�ĳ���6
   tccpara.A5=str2num(get(handles.edit_A5_1,'String'))*(10^6); 
   tccpara.phiA5=-str2num(get(handles.edit_A5_2,'String'));
   
   %A5�ĵ�λ��mm��nm��
   tccpara.C5=str2num(get(handles.edit_C5,'String'))*(10^6);
   
   %S5�ĵ�λ��mm��nm�� %�����ֵ��hrtem���滹û�г��֣���˽Ƕȹ�ϵ��ʱû�й�
   tccpara.S5=str2num(get(handles.edit_S5_1,'String'))*(10^6); tccpara.phiS5=-str2num(get(handles.edit_S5_2,'String'));
  
   %D5�ĵ�λ��mm��nm��%�����ֵ��hrtem���滹û�г��֣���˽Ƕȹ�ϵ��ʱû�й�
   tccpara.D5=str2num(get(handles.edit_D5_1,'String'))*(10^6); tccpara.phiD5=-str2num(get(handles.edit_D5_2,'String'));

   
   
   tccpara.A1x=tccpara.A1*cos(tccpara.phiA1/180*pi);   tccpara.A1y=tccpara.A1*sin(tccpara.phiA1/180*pi);
   %C1=focus;  focus���������defocs������ƫ���
   tccpara.A2x=tccpara.A2*cos(tccpara.phiA2/180*pi);   tccpara.A2y=tccpara.A2*sin(tccpara.phiA2/180*pi);
   tccpara.B2x=tccpara.B2*cos(tccpara.phiB2/180*pi);     tccpara.B2y=tccpara.B2*sin(tccpara.phiB2/180*pi);
   tccpara.A3x=tccpara.A3*cos(tccpara.phiA3/180*pi);   tccpara.A3y=tccpara.A3*sin(tccpara.phiA3/180*pi);
   tccpara.S3x=tccpara.S3*cos(tccpara.phiS3/180*pi);   tccpara.S3y=tccpara.S3*sin(tccpara.phiS3/180*pi);
   tccpara.C3=tccpara.Cs;
   tccpara.A4x=tccpara.A4*cos(tccpara.phiA4/180*pi);   tccpara.A4y=tccpara.A4*sin(tccpara.phiA4/180*pi);
   tccpara.D4x=tccpara.D4*cos(tccpara.phiD4/180*pi);   tccpara.D4y=tccpara.D4*sin(tccpara.phiD4/180*pi);
   tccpara.B4x=tccpara.B4*cos(tccpara.phiB4/180*pi);     tccpara.B4y=tccpara.B4*sin(tccpara.phiB4/180*pi);
   tccpara.A5x=tccpara.A5*cos(tccpara.phiA5/180*pi);   tccpara.A5y=tccpara.A5*sin(tccpara.phiA5/180*pi);
   tccpara.S5x=tccpara.S5*cos(tccpara.phiS5/180*pi);   tccpara.S5y=tccpara.S5*sin(tccpara.phiS5/180*pi);
   tccpara.C5=tccpara.C5;
   tccpara.D5x=tccpara.D5*cos(tccpara.phiD5/180*pi);   tccpara.D5y=tccpara.D5*sin(tccpara.phiD5/180*pi);

pp=1;
   
    

function x=myaperture(wave,sx,sy,kvector1,kvector2,shiftx,shifty,flag);%g_w is the radius  NOT dia
%sx��syΪ���������׿ռ䵥λ��ʸ������** ��-1��kvector1Ϊ�����Ĳ�ʸ���ȵ�Ƶ��kvector2Ϊ��ʸ��Ƶ������10 ��-1��shiftxΪƽ�ƵĲ�ʸ����������Ĵ��ƶ�һ���,flag��ʾ������ʽ��ȡkvector1��kvector2֮�����Ϣ�����ǲ�ȡ��֮�����Ϣ��
%�����Ƶ�˲���myaperture(newdd,1,1,0,kvector,0,0,1);
%����ȡ��ͨ��myaperture(newdd,1,1,kvector1,kvector,0,0,1);
%����ȡ��Ƶ��myaperture(newdd,1,1,0,kvector,0,0,0);
tx=shiftx/sx;
ty=shifty/sy;
 
[m,n]=size(wave);
uu=-round((n-1)/2):round(n/2)-1;
vv=-round((m-1)/2):round(m/2)-1;
[uu,vv]=meshgrid((uu-ty).*sy,(vv-tx).*sx);
uuvv=sqrt(uu.*uu+vv.*vv);
clear uu
clear vv
if flag==1
   wave(find(uuvv<kvector1 | uuvv>kvector2))=0;    %  �˲�
else
   wave(find(uuvv>=kvector1 & uuvv<=kvector2))=0;    %  �˲�
end
x=wave;
return;

function flag = getparaflag(handles)
%�õ����������ѡ�񣬷ֱ���p��l��n����
num = get(handles.popupmenu5,'value');
if num == 1
    flag = 'p';
elseif num == 2;
    flag = 'n';
elseif num ==3;
    flag = 'l';
end

% --- Executes on button press in pushbutton3.
function pushbutton3_Callback(hObject, eventdata, handles)
load allresul;
eachthick = allresul.eachthick;
handles.eachthick= eachthick;
slicethick = allresul.slicethick;

%�����м�������
display(strcat('Thickness of eachslice is ', num2str(eachthick),' angstrom, there are totally ', num2str(length(slicethick)), ' slices'));

handles.mid_slice_num = sort( str2num(get(handles.edit105,'string')));  % the mid slice's output
handles.mid_slice_num(find(handles.mid_slice_num>=length(slicethick))) = [];
display(strcat( 'Output the wave at the middle slices:', num2str(handles.mid_slice_num)))

paraflag = getparaflag(handles); 
%set(handles.pushbutton3,'enable','off');
%if paraflag == 'p'|| paraflag == 'l'  %����ʦ��ʽ 20210223
    [all_nuclear, all_nuc_ion, absorp_n, absorp_n_i, corr_peng_nuc, corr_peng_nuc_ion]=CommonPart_forsimulation_lobato_peng(hObject, eventdata, handles,paraflag);  %step0-step3�ƶ�������
%elseif paraflag == 'n' %����ʦ��ʽ������
     %�����peng����������õ�corr��ص��������󡣼���˼·���ǰ���������������Ƴ���֮����Ϊԭ�Ƴ��ĳ�ֵ
%end
    
%���ѭ��������ԭ�ӵ���������仯30��
%�ڲ�ѭ�����Ƿ��Ƿֲ�ЧӦ������ǣ�����Ҫ���㣬�ټ���ÿ��Ĺ��ס�

%all_nuclear ��һ����Ԫ�����࣬��2-4�������꣬��5����B���ӣ���6����ռ����
%step 4
%��ȡһЩ���������ѹ���������޹�
h_ba=(6.6255916e-34)/2/pi;
h=6.6255916e-34;
e=1.602102e-19;
me=9.1093807e-31;
c=2.9979251e+8;
vol=str2num(get(handles.edit_Vol,'string'));
handles.lambda=1.0e+9*h.*c./sqrt(e*vol*1000*(2*me*c*c+e*vol*1000));  %���㲨������λnm
handles.vol=vol*1000;

%step 5
%�޸�2020.04.29
if get(handles.radiobutton11,'value') ||  get(handles.radiobutton10,'value') ||  get(handles.radiobutton13,'value')  %stem�� or idpc or CBED�� 
   probesx=str2num(get(handles.edit4,'String'));
   [para_part1, U, V]=CommonPara_TEM(handles, probesx, probesx);  %��CTF��PhasePlate�ù̶��Ĵ�С����  
   para_part2 = readtccfromhandles_newTEM(hObject, handles, para_part1.lambda);   %ֻ��һ���������ѿ���е����ݶ���ȡ�ˣ�������֤�����ط���Ҫ��ȡ����ʱ����ͳһ�Ĵ���
    %������Ϊ�������֣�lambda��Ҫ����
   residual_aberr=WholeTCC_residual_newTEM(para_part2, U,V);  %����phase plate visualing the residual aberrations��kai'
    %-------------------------------------
   mytcc=WholeTCC2D_newTEM(para_part1, para_part2, U,V,2);
   
   myap=ones(probesx);  %�������������kmax�����涨��
   myap=myaperture(myap,1/(para_part1.sampling.*probesx),1/(para_part1.sampling.*probesx),0,para_part1.gmax*0.01/para_part1.lambda,0,0,1); %���10mrad���������1���ٳ���0.01��������rad���ٳ���nm��λ�Ĳ���
   num=length(mytcc(1,1,:));  %�����ж��ٸ�probe��Ҫɨ��ͼ��
   for i=1:num  %�������м䣬mytcc�������Ʊ����ˡ�
        guiyihua=sum(sum(abs(ifft2(ifftshift(myap.*mytcc(:,:,i)))).^2));  %20201226 add one aper to normalize 
        guiyihua=sqrt(guiyihua);  %20201231��Ҫ�󿪸��ţ�������probe����3.59��3.68ʽ
        handles.probe(:,:,i) =  fftshift(ifft2(ifftshift(myap.*mytcc(:,:,i))))./guiyihua;
   end
   handles.gfsf=para_part1.gfsf;
end
%__________�޸�2020.04.29 ��Ҫɾ��ԭ�е�STEM��׺����Щ����
if get(handles.radiobutton9,'value')  %HRTEM����  �������
    sy=handles.green_Nrow;sx=handles.green_Ncol;   %��ɫ���ĳߴ�
    [para_part1, U, V]=CommonPara_TEM(handles, sx, sy);  %��CTF��PhasePlate�ù̶��Ĵ�С����
    U=U-para_part1.tiltx;  %add 20210119 the incident wave is tilted but not the CTF is tilted
    V=V-para_part1.tilty;  %the tilted CTF is shown for CTF display. but in simulation, only a tilted incident beam is required.
 
    para_part2 = readtccfromhandles_newTEM(hObject, handles, para_part1.lambda);   %ֻ��һ���������ѿ���е����ݶ���ȡ�ˣ�������֤�����ط���Ҫ��ȡ����ʱ����ͳһ�Ĵ���
    %������Ϊ�������֣�lambda��Ҫ����
    residual_aberr=WholeTCC_residual_newTEM(para_part2, U,V);  %����phase plate visualing the residual aberrations
    %-------------------------------------
    mytcc=WholeTCC2D_newTEM(para_part1, para_part2, U,V,1);

    myap=ones(size(U));
    myap=myaperture(myap,1/(para_part1.sampling.*sx),1/(para_part1.sampling.*sy),0,para_part1.gmax,0,0,1);  
    
    handles.gfsf=para_part1.gfsf;
end


%step 6  ����
%�������ԭ��ʵ�ʵ�vibration����Ҫ���ɺܶ������
if get(handles.radiobutton2, 'value')   %if vibration is required
    allvib = str2num(get(handles.edit99,'string'));
else
    allvib = 1;
end

for vib = 1 : allvib;  %vibration
    %�������ԭ�ӵ�����ƫ��,��Ҫ��ÿ��ԭ�ӵ�Debye�йأ�֮����ӵ�ԭ�ӵ�������
    all_nuclear_copy=all_nuclear;
    all_nuc_ion_copy=all_nuc_ion;
    if paraflag=='n'
        corr_nuclear_copy = corr_peng_nuc;
        corr_peng_nuc_ion_copy = corr_peng_nuc_ion;
    end
        
    if ~isempty(all_nuclear) & get(handles.radiobutton2, 'value') % if vibration is required
        % make random positions
        % satisfy the condition that the atoms before special slices will
        % not be shifted to its slices
        % input atoms' parameters, thickness of each slice, the thickness
        % of the top, the slice number of middle slice.
        all_nuclear_copy = randshiftpos(all_nuclear_copy, eachthick, slicethick, handles.mid_slice_num); 
        if paraflag == 'n'
            if ~isempty(corr_peng_nuc)  %�����Ҫ�����Ƴ������������ֵĲ�����ԭ��λ�ø���һ��
                corr_nuclear_copy = corr_peng_nuc;
                corr_nuclear_copy(:,2:4) = all_nuclear_copy(:,2:4);
            end
        end
    end
    if ~isempty(all_nuc_ion) & get(handles.radiobutton2, 'value')
        all_nuc_ion_copy = randshiftpos(all_nuc_ion_copy, eachthick, slicethick, handles.mid_slice_num);
        if paraflag == 'n'
            if ~isempty(corr_peng_nuc_ion) %�����Ҫ�����Ƴ������������ֵĲ�����ԭ��λ�ø���һ��
                corr_peng_nuc_ion_copy = corr_peng_nuc_ion;
                corr_peng_nuc_ion_copy(:,2:4) = all_nuc_ion_copy(:,2:4);
            end
        end
    end
    
    %��Ҫ���߲���Ҫ����߶ȴ����Ĳ��,���һ���Ҫ���û�����ڷ�Χ�ڵ�ԭ�ӡ�    ~~~~~~~~~~~~~~~~~~~~~~~
    %improved multi-slice method
    flag=get(handles.checkbox1, 'value');  %�Ƿ���Ҫ����ԭ�ӻ����ڲ�ͬ�Ĳ�
    if flag==1
        multiceng=str2num(get(handles.edit98, 'string'));
        flag=multiceng;
        if paraflag == 'l'
            disp('Lobato parameters cannot be sliced into multiple slices in this program')
            flag=0;
        end
    end
    %�������ֲ�����ѡ��peng��lobato��
%     [ele_n, ele_n_i, absorp_n, absorp_n_i, series_n, series_n_i] = ...
%              CalAllEquilent20201106(all_nuclear_copy, all_nuc_ion_copy, absorp_n, absorp_n_i, eachthick, slicethick, flag);  
    DBmode=double(get(handles.radiobutton4,'value')|get(handles.radiobutton14,'value'));  %�κ�һ������1������1
    [ele_n, ele_n_i, absorp_n, absorp_n_i, series_n, series_n_i] = ...
             CalAllEquilent20210223_lobato_peng(all_nuclear_copy, all_nuc_ion_copy, absorp_n, absorp_n_i, eachthick, slicethick, flag, paraflag, DBmode);  
    
    ele_n_corr=[]; ele_n_i_corr=[]; corr_info=[]; series_n_corr=[]; series_n_i_corr=[];
    if paraflag == 'n'
        [ele_n_corr, ele_n_i_corr, corr_info, series_n_corr, series_n_i_corr] =  ...
            CalAllEquilent20210305_peng_corr(corr_nuclear_copy, corr_peng_nuc_ion_copy, eachthick, slicethick, flag);   
    end   
    
    %�õ���vib��ͼ��  ~~~~~~~~~~~~~~~~~~~~
    if vib==1
         if get(handles.radiobutton11,'value') || get(handles.radiobutton13,'value')  %stem��
             [myresul, midresul] = STEMsimulation(handles, ele_n, ele_n_i, absorp_n, absorp_n_i, series_n, series_n_i, ele_n_corr, ele_n_i_corr, corr_info, ... %��������peng��ϵ��������������Ƴ�����
    series_n_corr, series_n_i_corr);
         end
         if get(handles.radiobutton10,'value')  %cbed��, 
             myresul=CBEDsimulation(handles, ele_n, ele_n_i, absorp_n, absorp_n_i, series_n, series_n_i, ele_n_corr, ele_n_i_corr, corr_info, ... %��������peng��ϵ��������������Ƴ�����
    series_n_corr, series_n_i_corr);
         end
         if get(handles.radiobutton9,'value')  %hrtem��
             myresul=HRTEMsimulation(handles, ele_n, ele_n_i, absorp_n, absorp_n_i, series_n, series_n_i, mytcc, ele_n_corr, ele_n_i_corr, corr_info, ... %��������peng��ϵ��������������Ƴ�����
    series_n_corr, series_n_i_corr);
             
         end
     else
         vib
         if get(handles.radiobutton11,'value') || get(handles.radiobutton13,'value')  %stem��
             [myresultemp, midresultemp] = STEMsimulation(handles, ele_n, ele_n_i, absorp_n, absorp_n_i, series_n, series_n_i, ele_n_corr, ele_n_i_corr, corr_info, ... %��������peng��ϵ��������������Ƴ�����
    series_n_corr, series_n_i_corr);
             myresul=myresul+myresultemp;
             midresul = midresul+midresultemp; 
         end
         if get(handles.radiobutton10,'value')  %cbed��
             myresul=myresul+CBEDsimulation(handles, ele_n, ele_n_i, absorp_n, absorp_n_i, series_n, series_n_i,ele_n_corr, ele_n_i_corr, corr_info, ... %��������peng��ϵ��������������Ƴ�����
    series_n_corr, series_n_i_corr);
         end
         if get(handles.radiobutton9,'value')  %hrtem��
             myresul=myresul+HRTEMsimulation(handles, ele_n, ele_n_i, absorp_n, absorp_n_i, series_n, series_n_i, mytcc, ele_n_corr, ele_n_i_corr, corr_info,series_n_corr, series_n_i_corr);
         end
      end
end

disp(strcat('Results will be saved in the folder:', handles.saveresult));
nname = get(handles.edit106,'string');

allstemnum = 1; %����stemͼ��ĸ���
steminfo = [];
if get(handles.radiobutton11,'value')  %stem�� ��Ҫ����reshape�����С
    aper_dect = str2num(get(handles.edit79,'string'));%������չ��� 

   %�м��� & �����
   if isempty(midresul)
       midresul=myresul(:);
   else
      midresul(:,end+1)=myresul(:);
   end
   handles.mid_slice_num = [handles.mid_slice_num, length(slicethick)];
   midimg_result_temp = zeros( handles.width_red , handles.hight_red, length(myresul(1,:)));
   for i=1:length(handles.mid_slice_num)
       midimg_result_temp(:) = midresul(:,i)/probesx/probesx/allvib;
              
       for j=1:length(myresul(1,:))
            figure;imshow(midimg_result_temp(:,:,j).', 'XData',...
           str2num(get(handles.edit17,'string')) + ([1:str2num(get(handles.edit20,'string'))]-1).*str2num(get(handles.edit77,'string')).*str2num(get(handles.edit75,'string')), ...
           'YData', str2num(get(handles.edit18,'string')) + ([1:str2num(get(handles.edit19,'string'))]-1).*str2num(get(handles.edit77,'string')).*str2num(get(handles.edit75,'string')),...
          'DisplayRange',[]);axis on;xlabel('Angstrom');
      
            distitle = strcat('STEM images (thickness:', ...
                num2str(handles.mid_slice_num(i)*eachthick) , ...
                '; Dectector: ', ...
                num2str(aper_dect(j,1)), ...
                'to ', num2str(aper_dect(j,2))   );colorbar
            title( distitle )
            
            discontent = strcat('Mid result _W', num2str(handles.width_red),' * H', num2str(handles.hight_red), '_Float Data in',...
                strcat(handles.saveresult, '\',nname ,'_', 'STEM_', num2str(handles.mid_slice_num(i)*eachthick), ...
                'slice_Dect_',num2str(round(aper_dect(j,1))),'_',num2str(round(aper_dect(j,2))),'.dat'))
            disp(discontent);
            
            fname = strcat(handles.saveresult, '\',nname ,'_', 'STEM_', num2str(handles.mid_slice_num(i)*eachthick), ...
                'slice_Dect_',num2str(round(aper_dect(j,1))),'_',num2str(round(aper_dect(j,2))),'.dat')
            fid = fopen(fname,'w');
            fwrite(fid, midimg_result_temp(:,:,j), 'float'); 
            fclose(fid);
            
            steminfo(allstemnum).fname = fname;
            steminfo(allstemnum).discontent = discontent;
            steminfo(allstemnum).title = distitle;
            allstemnum = allstemnum+1;
       end
   end
   
end

if get(handles.radiobutton9,'value')  %hrtem�� ��Ҫ����reshape�����С
    for i=1:length(myresul(1,1,:))
      figure;imshow(myresul(:,:,i)./allvib, 'XData',...
           str2num(get(handles.edit17,'string')) + ([1:str2num(get(handles.edit20,'string'))]-1).*str2num(get(handles.edit77,'string')).*str2num(get(handles.edit75,'string')), ...
           'YData', str2num(get(handles.edit18,'string')) + ([1:str2num(get(handles.edit19,'string'))]-1).*str2num(get(handles.edit77,'string')).*str2num(get(handles.edit75,'string')),...
          'DisplayRange',[]);axis on;xlabel('Angstrom');colorbar
      
      disp(strcat('Result _W', num2str(handles.width_red),' * H', num2str(handles.hight_red), ...
          '_Float Data in', strcat(handles.saveresult, '\', nname ,'_', 'HRTEM','.dat')));
      fid = fopen(strcat(handles.saveresult, '\', nname ,'_', 'HRTEM','.dat'),'w');
      fwrite(fid, myresul(:,:,i).'/allvib, 'float'); 
      fclose(fid);
    end
end
if get(handles.radiobutton10,'value')  %CBED�� ��Ҫ����reshape�����С
    for i=1:length(myresul(1,1,:))
      figure;imshow(myresul(:,:,i)./allvib,[])
      
      disp(strcat('Result _W', num2str(handles.CBEDprobesx),' * H', num2str(handles.CBEDprobesx), ...
          '_Float Data in', strcat(handles.saveresult, '\', nname ,'_', 'CBED','.dat')));
      fid = fopen(strcat(handles.saveresult, '\', nname ,'_', 'CBED','.dat'),'w');
      fwrite(fid, myresul(:,:,i).'/allvib, 'float'); 
      fclose(fid);
    end
end

% step=2;
% sampling=0.1;
% [x,y] = meshgrid( ([1:newy]-mycen(newy))/(newy*sampling*step/handles.conv_sampling) ,([1:newx]-mycen(newx))/(newx*sampling*step/handles.conv_sampling));
%         ss = 0.7/(sqrt(-log(0.5))*2);
%Kepeng
if get(handles.radiobutton13,'value')  %idpc��
    aper_dect = str2num(get(handles.edit79,'string'));%������չ��� 

    [kx,ky]=meshgrid((-handles.width_red/2:handles.width_red/2-1)./(handles.width_red.*para_part1.sampling*str2num(get(handles.edit77,'string'))),(-handles.hight_red/2:handles.hight_red/2-1)/(handles.hight_red.*para_part1.sampling*str2num(get(handles.edit77,'string'))));
    k2=kx.^2+ky.^2;  %��ֵ��Χ��һ����Ҫ���
    ss = 0.1/(sqrt(-log(0.5))*2);

   %�м��� & �����
     if isempty(midresul)
         midresul=myresul(:);
     else
         midresul(:,end+1)=myresul(:);
     end
     handles.mid_slice_num = [handles.mid_slice_num, length(slicethick)];
     midimg_result_temp = zeros( handles.width_red,handles.hight_red, length(myresul(1,:)));
     for i=1:length(handles.mid_slice_num)
          midimg_result_temp(:) = midresul(:,i);
          
          ta=zeros(handles.width_red,handles.hight_red);    
          for j=1:(length(myresul(1,:))/4)
              ta(:)=midimg_result_temp(:,:,4*(j-1)+1)/probesx/probesx/allvib;
              aa(:,:,j)=ta';
              ta(:)=midimg_result_temp(:,:,4*(j-1)+2)/probesx/probesx/allvib;
              bb(:,:,j)=ta';
              ta(:)=midimg_result_temp(:,:,4*(j-1)+3)/probesx/probesx/allvib;
              cc(:,:,j)=ta';
              ta(:)=midimg_result_temp(:,:,4*(j-1)+4)/probesx/probesx/allvib;
              dd(:,:,j)=ta';
%               ta(:)=midimg_result_temp(:,:,4*(j-1)+1)/probesx/probesx/allvib;
%               t=fftshift(fft2(ta.'));aa(:,:,j)=ifft2(ifftshift(t.*exp(-((pi*ss)^2*k2))));
%               ta(:)=midimg_result_temp(:,:,4*(j-1)+2)/probesx/probesx/allvib;
%               t=fftshift(fft2(ta.'));bb(:,:,j)=ifft2(ifftshift(t.*exp(-((pi*ss)^2*k2))));
%               ta(:)=midimg_result_temp(:,:,4*(j-1)+3)/probesx/probesx/allvib;
%               t=fftshift(fft2(ta.'));cc(:,:,j)=ifft2(ifftshift(t.*exp(-((pi*ss)^2*k2))));
%               ta(:)=midimg_result_temp(:,:,4*(j-1)+4)/probesx/probesx/allvib;
%               t=fftshift(fft2(ta.'));dd(:,:,j)=ifft2(ifftshift(t.*exp(-((pi*ss)^2*k2))));
       
        
%                 px=(fftshift(fft2(bb(:,:,j)-dd(:,:,j))).*kx)*(2*pi*sqrt(-1)); 
%                 px(handles.hight_red/2+1,handles.width_red/2+1)=mean(mean(bb(:,:,j)-dd(:,:,j)));
%                 py=(fftshift(fft2(cc(:,:,j)-aa(:,:,j))).*ky)*(2*pi*sqrt(-1)); 
%                 py(handles.hight_red/2+1,handles.width_red/2+1)=mean(mean(cc(:,:,j)-aa(:,:,j)));
%                 px=real( ifft2(ifftshift(px)) );
%                 py=real( ifft2(ifftshift(py)) );
                px = bb(:,:,j)-dd(:,:,j);
                py = cc(:,:,j)-aa(:,:,j);
                
                fid = fopen(strcat(handles.saveresult, '\',nname ,'_', 'aa_',num2str(handles.mid_slice_num(i)),'slice_Dect_', num2str(round(aper_dect(j,1))), ...
                    '_', num2str(round(aper_dect(j,2))),'.dat'),'w');
                fwrite(fid, aa.', 'float'); 
                fclose(fid);
                fid = fopen(strcat(handles.saveresult, '\',nname ,'_', 'bb_',num2str(handles.mid_slice_num(i)),'slice_Dect_', num2str(round(aper_dect(j,1))), ...
                    '_', num2str(round(aper_dect(j,2))),'.dat'),'w');
                fwrite(fid, bb.', 'float'); 
                fclose(fid);
                fid = fopen(strcat(handles.saveresult, '\',nname ,'_', 'cc_',num2str(handles.mid_slice_num(i)),'slice_Dect_', num2str(round(aper_dect(j,1))), ...
                    '_', num2str(round(aper_dect(j,2))),'.dat'),'w');
                fwrite(fid, cc.', 'float'); 
                fclose(fid);
                fid = fopen(strcat(handles.saveresult, '\',nname ,'_', 'dd_',num2str(handles.mid_slice_num(i)),'slice_Dect_', num2str(round(aper_dect(j,1))), ...
                    '_', num2str(round(aper_dect(j,2))),'.dat'),'w');
                fwrite(fid, dd.', 'float'); 
                fclose(fid);
                
                
                
                figure;subplot(2,2,1);imshow(px, 'XData',...
                      str2num(get(handles.edit17,'string')) + ([1:str2num(get(handles.edit20,'string'))]-1).*str2num(get(handles.edit77,'string')).*str2num(get(handles.edit75,'string')), ...
                      'YData', str2num(get(handles.edit18,'string')) + ([1:str2num(get(handles.edit19,'string'))]-1).*str2num(get(handles.edit77,'string')).*str2num(get(handles.edit75,'string')),...
                      'DisplayRange',[]);axis on;xlabel('Angstrom');
                title('x vector');subplot(2,2,2);imshow(py, 'XData',...
                      str2num(get(handles.edit17,'string')) + ([1:str2num(get(handles.edit20,'string'))]-1).*str2num(get(handles.edit77,'string')).*str2num(get(handles.edit75,'string')), ...
                      'YData', str2num(get(handles.edit18,'string')) + ([1:str2num(get(handles.edit19,'string'))]-1).*str2num(get(handles.edit77,'string')).*str2num(get(handles.edit75,'string')),...
                      'DisplayRange',[]);axis on;xlabel('Angstrom');
                title('y vector');
                [sx,sy]=size(px);newp=zeros(sx,sy,3);newp(:,:,1)=px;newp(:,:,2)=py;
                subplot(2,2,3);image(newp*255);axis on;xlabel('Angstrom');
                subplot(2,2,4);imshow(-(px).^2+(py).^2, 'XData',...
                       str2num(get(handles.edit17,'string')) + ([1:str2num(get(handles.edit20,'string'))]-1).*str2num(get(handles.edit77,'string')).*str2num(get(handles.edit75,'string')), ...
                       'YData', str2num(get(handles.edit18,'string')) + ([1:str2num(get(handles.edit19,'string'))]-1).*str2num(get(handles.edit77,'string')).*str2num(get(handles.edit75,'string')),...
                       'DisplayRange',[]);axis on;xlabel('Angstrom');
                title('x^2+y^2');    
                title( strcat('DDPC (thick:', ...
                      num2str(handles.mid_slice_num(i)*eachthick) , ...
                      '; Dectector: ', ...
                      num2str(aper_dect(j,1)), ...
                     'to ', num2str(aper_dect(j,2))   ))
                disp(strcat('Result _W', num2str(handles.width_red),' * H', num2str(handles.hight_red),...
                    '_Float Data in', strcat(handles.saveresult, '\',nname ,'_', 'DDPC_', num2str(handles.mid_slice_num(i)),'slice_Dect_',num2str(round(aper_dect(j,1))), ...
                    '_', num2str(round(aper_dect(j,2))),'.dat')));
                fid = fopen(strcat(handles.saveresult, '\',nname ,'_', 'DDPCx_',num2str(handles.mid_slice_num(i)),'slice_Dect_', num2str(round(aper_dect(j,1))), ...
                    '_', num2str(round(aper_dect(j,2))),'.dat'),'w');
                fwrite(fid, px.', 'float'); 
                fclose(fid);
                fid = fopen(strcat(handles.saveresult, '\',nname ,'_', 'DDPCy_', num2str(handles.mid_slice_num(i)),'slice_Dect_',num2str(round(aper_dect(j,1))), ...
                     '_', num2str(round(aper_dect(j,2))),'.dat'),'w');
                fwrite(fid, py.', 'float'); 
                fclose(fid);
                 
                 
                p=(fftshift(fft2(px)).*kx+fftshift(fft2(py)).*ky)./(2*pi*sqrt(-1)*k2); 
                p(handles.hight_red/2+1,handles.width_red/2+1)=(mean(mean(px))+mean(mean(py)));
                p=real( ifft2(ifftshift(p)) );
                
                figure;imshow(p, 'XData',...
                      str2num(get(handles.edit17,'string')) + ([1:str2num(get(handles.edit20,'string'))]-1).*str2num(get(handles.edit77,'string')).*str2num(get(handles.edit75,'string')), ...
                      'YData', str2num(get(handles.edit18,'string')) + ([1:str2num(get(handles.edit19,'string'))]-1).*str2num(get(handles.edit77,'string')).*str2num(get(handles.edit75,'string')),...
                      'DisplayRange',[]);axis on;xlabel('Angstrom');colorbar;
                title( strcat('IDPC (thick:', ...
                      num2str(handles.mid_slice_num(i)*eachthick) , ...
                      '; Dectector: ', ...
                      num2str(aper_dect(j,1)), ...
                     'to ', num2str(aper_dect(j,2))   ))
                disp(strcat('Result _W', num2str(handles.width_red),' * H', num2str(handles.hight_red), ...
                     '_Float Data in', strcat(handles.saveresult, '\',nname ,'_', 'IDPC_', num2str(handles.mid_slice_num(i)*eachthick), 'slice_Dect_',num2str(round(aper_dect(j,1))), ...
                     '_', num2str(round(aper_dect(j,2))),'.dat')));
                fid = fopen(strcat(handles.saveresult, '\',nname ,'_', 'IDPC_', num2str(handles.mid_slice_num(i)*eachthick), 'slice_Dect_',num2str(round(aper_dect(j,1))), ...
                     '_', num2str(round(aper_dect(j,2))),'.dat'),'w');
                fwrite(fid, p.', 'float'); 
                fclose(fid);
            
%                 figure;imshow((aa(:,:,j)+cc(:,:,j)+bb(:,:,j)+dd(:,:,j)), 'XData',...
%                        str2num(get(handles.edit17,'string')) + ([1:str2num(get(handles.edit20,'string'))]-1).*str2num(get(handles.edit77,'string')).*str2num(get(handles.edit75,'string')), ...
%                        'YData', str2num(get(handles.edit18,'string')) + ([1:str2num(get(handles.edit19,'string'))]-1).*str2num(get(handles.edit77,'string')).*str2num(get(handles.edit75,'string')),...
%                        'DisplayRange',[]);axis on;xlabel('Angstrom');colorbar;
%               
%                 title( strcat('Traditional STEM (thick:', ...
%                       num2str(handles.mid_slice_num(i)*eachthick) , ...
%                       '; Dectector: ', ...
%                       num2str(aper_dect(j,1)), ...
%                      'to ', num2str(aper_dect(j,2))   ))
%                  disp(strcat('Result _W', num2str(handles.width_red),' * H', num2str(handles.hight_red), ...
%                       '_Float Data in', strcat(handles.saveresult, '\',nname ,'_', 'IDPC_STEM', num2str(handles.mid_slice_num(i)*eachthick), 'slice_Dect_', num2str(round(aper_dect(j,1))), ...
%                       '_', num2str(round(aper_dect(j,2))),'.dat')));
%                  fid = fopen(strcat(handles.saveresult, '\',nname ,'_', 'IDPC_STEM_', num2str(handles.mid_slice_num(i)*eachthick), 'slice_Dect_',num2str(round(aper_dect(j,1))), ...
%                         '_', num2str(round(aper_dect(j,2))),'.dat'),'w');
%                  fwrite(fid, aa(:,:,j)+cc(:,:,j)+bb(:,:,j)+dd(:,:,j), 'float'); 
%                  fclose(fid);
          end
   end

end
    
 
%Kepeng
% handles.steminfo = steminfo; %save all file information to do convolution
% guidata(hObject, handles);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%KepengWu1124   begin
function APER =makeAPER(handles);
aper_dect = str2num(get(handles.edit79,'string'));%������չ���
[px, py]=meshgrid( (-handles.probesx/2:(handles.probesx/2-1))./(handles.probesx*handles.sampling), ...
    (-handles.probesx/2:(handles.probesx/2-1))./(handles.probesx*handles.sampling));
p2_probe=(px.^2+py.^2);
aper_range = aper_dect*0.001/(handles.lambda*10); %1/A ��λ
num=length(aper_dect(:,1));
APER=zeros(handles.probesx, handles.probesx,num);

[theta, rho] =  cart2pol(px,py);  % range of theta is from -pi to pi
for i=1:num  %�������չ���
    temp=zeros(handles.probesx, handles.probesx);
    temp(find(p2_probe>=aper_range(i,1).^2 & p2_probe<=aper_range(i,2).^2 )) = 1;

    if aper_dect(i,3) >180  %��0-360����д��Χת�䵽-180 ��180
        aper_dect(i,3) = aper_dect(i,3)-180; 
    end
    if aper_dect(i,3) >= -180 & aper_dect(i,3) < -90
        dectvalue = aper_dect(i,3);
    elseif aper_dect(i,3) >= -90 & aper_dect(i,3) < 0
        dectvalue = aper_dect(i,3)-90;
    elseif aper_dect(i,3) >= 0 & aper_dect(i,3) < 90
        dectvalue = aper_dect(i,3)-180;
    elseif aper_dect(i,3) >= 90 & aper_dect(i,3) < 180
        dectvalue = aper_dect(i,3)-270;
    end
    tempaper = zeros(handles.probesx, handles.probesx);
    tempaper( find(theta>= dectvalue*pi/180 & theta<(dectvalue+90)*pi/180) ) = 1;
    APER(:,:,4*(i-1)+1) = tempaper;
    tempaper = zeros(handles.probesx, handles.probesx);
    tempaper( find(theta>= (dectvalue+90)*pi/180 & theta<(dectvalue+180)*pi/180) ) = 1;
    APER(:,:,4*(i-1)+2) = tempaper;
    tempaper = zeros(handles.probesx, handles.probesx);
    tempaper( find(theta>= (dectvalue+180)*pi/180 & theta<(dectvalue+270)*pi/180) ) = 1;
    APER(:,:,4*(i-1)+3) = tempaper;
    
    APER(:,:,4*(i-1)+4) = 1- APER(:,:,4*(i-1)+1) -APER(:,:,4*(i-1)+2) - APER(:,:,4*(i-1)+3);

    APER(:,:,4*(i-1)+1) = temp.* APER(:,:,4*(i-1)+1);  %���ϴ�ͨ��
    APER(:,:,4*(i-1)+2) = temp.* APER(:,:,4*(i-1)+2);
    APER(:,:,4*(i-1)+3) = temp.* APER(:,:,4*(i-1)+3);
    APER(:,:,4*(i-1)+4) = temp.* APER(:,:,4*(i-1)+4);
    % APER(:,:,4*(i-1)+1)=temp;
   
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  KepengWu1124 end
function myresul=HRTEMsimulation(handles, ele_n, ele_n_i, absorp_n, absorp_n_i, series_n, series_n_i, mytcc,...
    ele_n_corr, ele_n_i_corr, corr_info, ... %��������peng��ϵ��������������Ƴ�����
    series_n_corr, series_n_i_corr)
%��Ҫ����green�ĳߴ磬�Լ�����ԭ�ӵ����λ�ã�����������Լ�ƴ�ջ�ԭͼʱ�������
h_ba=(6.6255916e-34)/2/pi;
h=6.6255916e-34;
e=1.602102e-19;
me=9.1093807e-31;
c=2.9979251e+8;
PARAMETER=2*pi*(h_ba)*(h_ba)/(e*me*1.0E-10*1.0E-10)/(handles.green_Nrow*handles.green_Ncol*handles.sampling*handles.sampling);  %Angstrom��λ
sigma=2*pi/(handles.lambda*10*handles.vol/1000)*(me*c*c+e*handles.vol)/(2*me*c*c+e*handles.vol);  %��λ����һ����/kV*A

%APERTURE  % �������󣬱�֤�����������wrapЧӦ
APERTURE=ones(handles.green_Ncol, handles.green_Nrow);
[tx, ty]= meshgrid(-handles.green_Nrow/2:handles.green_Nrow/2-1, -handles.green_Ncol/2:handles.green_Ncol/2-1);
minn=min(handles.green_Nrow/2,handles.green_Ncol/2);
%APERTURE( find((tx./(handles.green_Nrow)).^2+(ty/(handles.green_Ncol)).^2>0.25/4) )=0;
%�ش�Ľ���20210225�������Ƶ�ľ�����Ϣ��ʧ��
clear tx ty

paraflag = getparaflag(handles);
if get(handles.GPURB,'value')  %ʹ��GPU������
    %GPU�����Ƿ���
    %%��ɫ�ĳߴ磬��ɫ����ʸ��gx���ڼ���ԭ��λ�õģ�s2ʹ��peng�Ĳ�������ԭ���Ƴ��ֲ�
    gncol = handles.green_Ncol; gnrow = handles.green_Nrow;
    green=max(handles.green_Nrow, handles.green_Ncol);
    handles.green_Ncol=green; handles.green_Nrow=green;
    [gx_green,gy_green]=meshgrid((-green/2:(green/2-1))./(green*handles.sampling), ...
        (-green/2:(green/2-1))./(green*handles.sampling)); %REPRO��λ��1/A)^2,����(1/nm)^2
    sx_green=gx_green/2;  % peng ��s����
    sy_green=gy_green/2;
    s2_green=sx_green.^2+sy_green.^2;
    APERTURE=ones(green, green);
    
    if paraflag=='n'
        disp('Calculating the correction for Peng''s scattering factor')
        corr_info_matrix=0.*s2_green;%����ֵ
        for i = 1:length(corr_info(:,1))  %��ÿ�ֵ������Ƴ����������������٣���������Ԫ�ظ����൱�ľ���ͺ�
            r2 = s2_green;
            r2(find(r2>=corr_info(i,5).^2 & r2<corr_info(i,6).^2))=(sin(pi * ((sqrt(r2(find(r2>=corr_info(i,5).^2 & r2<corr_info(i,6).^2)) )-corr_info(i,5))./(corr_info(i,6)-corr_info(i,5))-0.5) )+1)/2;
            r2(find(r2<corr_info(i,5).^2))=0;
            r2(find(r2>=corr_info(i,6).^2))=1;
            corr_info_matrix(:,:,i) = r2.*( corr_info(i,1).*exp(-s2_green.*corr_info(i,3)) + corr_info(i,2).*exp(-s2_green.*corr_info(i,4)));                           
        end  %gpu,����������������һЩ
    end
    
    handles.probe=[];
    APER=[];
    psf_fft=[];
    handles.probe_ingreenNrow=[]; handles.probe_ingreenNcol=[]; handles.probestep=[];
    tic
    cuda_stem_potential;
    toc
    
    handles.green_Ncol=gncol; handles.green_Nrow=gnrow;
    potential_temp = zeros(green,green, max(length(series_n), length(series_n_i)) );
    potential_temp(:) = potentialx(:) + sqrt(-1)*potentialy(:);
    
    for i=1:length(potential_temp(1,1,:)); %��Ҫת��
        potential_temp(:,:,i) = potential_temp(:,:,i).';
    end
    potential = potential_temp(1:handles.green_Ncol,1:handles.green_Nrow, :);
    clear potential_temp;
end
%CPU�����Ǿ���
    %%��ɫ�ĳߴ磬��ɫ����ʸ��gx���ڼ���ԭ��λ�õģ�s2ʹ��peng�Ĳ�������ԭ���Ƴ��ֲ�
[gx_green,gy_green]=meshgrid((-handles.green_Nrow/2:(handles.green_Nrow/2-1))./(handles.green_Nrow*handles.sampling), ...
        (-handles.green_Ncol/2:(handles.green_Ncol/2-1))./(handles.green_Ncol*handles.sampling)); %REPRO��λ��1/A)^2,����(1/nm)^2
    sx_green=gx_green/2;  % peng ��s����
    sy_green=gy_green/2;
    s2_green=sx_green.^2+sy_green.^2;
    

if get(handles.CPURB,'value')  %ʹ��CPU������ 
    
%����ÿ����Ƴ�
% potential=GetPotential4AllSlice_multicore_lobato_peng(handles.green_Ncol, handles.green_Nrow,... 
%     ele_n, absorp_n, ....�������ӻ���ԭ�ӵĵ��Ժ�����
%     ele_n_i, absorp_n_i, ... %ԭ�ӻ�ԭ��+���ӣ����Ի�������
%     series_n, series_n_i, ...  %ԭ�����д���   
%     s2_green, gx_green, gy_green, ...
%     sigma, PARAMETER, APERTURE, paraflag);   %Ϊ��STEM���㲻��������ֻ���뵽HRTEM��CBED��
APERTURE=ones(handles.green_Ncol, handles.green_Nrow);
potential=GetPotential4AllSlice_multicore_lotabo_peng_corr(handles.green_Ncol, handles.green_Nrow,... 
    ele_n, absorp_n, ....�������ӻ���ԭ�ӵĵ��Ժ�����  %��������������peng�Ĳ���
    ele_n_i, absorp_n_i, ... %ԭ�ӻ�ԭ��+���ӣ����Ի�������
    series_n, series_n_i, ...  %ԭ�����д���   
    ele_n_corr, ele_n_i_corr, corr_info, ... %�������peng��ϵ��������������Ƴ�����
    series_n_corr, series_n_i_corr,...
    s2_green, gx_green, gy_green, ...
    sigma, PARAMETER, APERTURE, paraflag);   %Ϊ��STEM���㲻��������ֻ���뵽HRTEM��CBED��
end


%�������
myap=ones(size(gx_green));
myap=myaperture(myap,1/(handles.sampling.*handles.green_Ncol),1/(handles.sampling.*handles.green_Nrow),0,str2num(get(handles.edit_Ape,'string'))*0.1,0,0,1);  %���ĵ�λ

%��ĳ�ʼ��
myinten=zeros(size(gx_green));

%add 20210119
Vol=str2num(get(handles.edit_Vol,'String')); %��ѹ
lambda=1.0e+9*h.*c./sqrt(e*Vol*1000*(2*me*c*c+e*Vol*1000));  %���㲨������λnm
tilt=str2num(get(handles.edit_Tilt_1,'String'))*10^(-3);
phitilt=-str2num(get(handles.edit_Tilt_2,'String'));  % ���ȵ�λ ���� 10mrad  
tiltx=tilt/lambda*cos(phitilt/180*pi); %���㵽���׸�ʽ�ĵ�λ 1/nm
tilty=tilt/lambda*sin(phitilt/180*pi);
tiltxnum=tiltx*handles.sampling*0.1*handles.green_Ncol;  %nm��sampling ���Ļ���
tiltynum=tilty*handles.sampling*0.1*handles.green_Ncol;

%����������  %20210119 ������б�Ĳ�����----------------

%   ������propogation;��zheight�й�
p2=((gx_green+tiltx*0.1).^2+(gy_green+tilty*0.1).^2);  %������
psf_fft=exp(-sqrt(-1)*pi*handles.lambda*10*handles.eachthick*p2);   %���ڵ�λA���㡣
psf_fft(find(p2>1/(16*handles.sampling*handles.sampling)))=0;  %��������
psf_fft=ifftshift(psf_fft);

myfftwave = zeros(size(gx_green)); 
if tilt > 0
     myfftwave(mycen(handles.green_Ncol)-round(tiltynum), mycen(handles.green_Nrow)-round(tiltxnum)) = 1;
else
     myfftwave(mycen(handles.green_Ncol), mycen(handles.green_Nrow)) = 1;  %��ֵΪƽ�沨
end
   %this wave shift is the same as the U*U+V*V before about line 885 U=U-para_part1.tiltx;   
mywave=ifft2(ifftshift(myfftwave))*handles.green_Nrow*handles.green_Ncol;  %���벨����
%----------end 20210119

for k=1:length(potential(1,1,:));
    if rem(k,10)==0
        pp=1;
    end
        mywave =  mywave.*potential(:,:,k);
        mywave=fft2(mywave);   %ע�⣬��ʡ��һ��fftshift
        mywave=mywave.*psf_fft;
        mywave=ifft2(mywave);
end
%�沨����-----add for Yao fenfa
[sx, sy]= size(gx_green);
disp(strcat('Wave saved in HRTEM_exitwave.dat (complex 8), H',num2str(sx), '*W',num2str(sy)))
a=zeros(2*sx,sy);
a(1:2:end,:)=real(mywave);
a(2:2:end,:)=imag(mywave);
fid = fopen('HRTEM_exitwave.dat','w');
fwrite(fid, a, 'float')
fclose(fid);

mywave=fftshift(fft2(mywave));

    
for nn=1:length(handles.gfsf);
    myinten=myinten+handles.gfsf(nn).*abs(ifft2(ifftshift(mywave.*mytcc(:,:,nn).*myap))).^2;
end
myresul=myinten(handles.HRTEM_ingreenNcol: handles.HRTEM_ingreenNcol+handles.hight_red-1, ...
                handles.HRTEM_ingreenNrow: handles.HRTEM_ingreenNrow+handles.width_red-1);
return;


function myresul=CBEDsimulation(handles, ele_n, ele_n_i, absorp_n, absorp_n_i, series_n, series_n_i, ele_n_corr, ele_n_i_corr, corr_info, ... %��������peng��ϵ��������������Ƴ�����
    series_n_corr, series_n_i_corr);
%��Ҫ����green�ĳߴ磬�Լ�����ԭ�ӵ����λ�ã�����������Լ�ƴ�ջ�ԭͼʱ�������
h_ba=(6.6255916e-34)/2/pi;
h=6.6255916e-34;
e=1.602102e-19;
me=9.1093807e-31;
c=2.9979251e+8;
PARAMETER=2*pi*(h_ba)*(h_ba)/(e*me*1.0E-10*1.0E-10)/(handles.green_Nrow*handles.green_Ncol*handles.sampling*handles.sampling);  %Angstrom��λ
sigma=2*pi/(handles.lambda*10*handles.vol/1000)*(me*c*c+e*handles.vol)/(2*me*c*c+e*handles.vol);  %��λ����һ����/kV*A

%APERTURE  % �������󣬱�֤�����������wrapЧӦ
APERTURE=ones(handles.green_Ncol, handles.green_Nrow);
[tx, ty]= meshgrid(-handles.green_Nrow/2:handles.green_Nrow/2-1, -handles.green_Ncol/2:handles.green_Ncol/2-1);
minn=min(handles.green_Nrow/2,handles.green_Ncol/2);
%APERTURE( find((tx./(handles.green_Nrow)).^2+(ty/(handles.green_Ncol)).^2>0.25/4) )=0;
clear tx ty

%%��ɫ�ĳߴ磬��ɫ����ʸ��gx���ڼ���ԭ��λ�õģ�s2ʹ��peng�Ĳ�������ԭ���Ƴ��ֲ�
[gx_green,gy_green]=meshgrid((-handles.green_Nrow/2:(handles.green_Nrow/2-1))./(handles.green_Nrow*handles.sampling), ...
    (-handles.green_Ncol/2:(handles.green_Ncol/2-1))./(handles.green_Ncol*handles.sampling)); %REPRO��λ��1/A)^2,����(1/nm)^2
sx_green=gx_green/2;  % peng ��s����
sy_green=gy_green/2;
s2_green=sx_green.^2+sy_green.^2;



%   ������propogation;��zheight�й�
[px, py]=meshgrid( (-handles.CBEDprobesx/2:(handles.CBEDprobesx/2-1))./(handles.CBEDprobesx*handles.sampling), ...
    (-handles.CBEDprobesx/2:(handles.CBEDprobesx/2-1))./(handles.CBEDprobesx*handles.sampling));
p2_probe=(px.^2+py.^2);  %������
psf_fft=exp(-sqrt(-1)*pi*handles.lambda*10*handles.eachthick*p2_probe);   %���ڵ�λA���㡣
%psf_fft(find(p2_probe>1/(16*handles.sampling*handles.sampling)))=0;  %��������
psf_fft(find(p2_probe>1/(9*handles.sampling*handles.sampling)))=0;  %����2/3�Ĺ���
psf_fft=ifftshift(psf_fft);

paraflag = getparaflag(handles);
if get(handles.GPURB,'value')  %ʹ��GPU������
    %GPU�����Ƿ���
    if paraflag=='n'
        disp('Calculating the correction for Peng''s scattering factor')
        corr_info_matrix=0.*s2_green;%����ֵ
        for i = 1:length(corr_info(:,1))  %��ÿ�ֵ������Ƴ����������������٣���������Ԫ�ظ����൱�ľ���ͺ�
            r2 = s2_green;
            r2(find(r2>=corr_info(i,5).^2 & r2<corr_info(i,6).^2))=(sin(pi * ((sqrt(r2(find(r2>=corr_info(i,5).^2 & r2<corr_info(i,6).^2)) )-corr_info(i,5))./(corr_info(i,6)-corr_info(i,5))-0.5) )+1)/2;
            r2(find(r2<corr_info(i,5).^2))=0;
            r2(find(r2>=corr_info(i,6).^2))=1;
            corr_info_matrix(:,:,i) = r2.*( corr_info(i,1).*exp(-s2_green.*corr_info(i,3)) + corr_info(i,2).*exp(-s2_green.*corr_info(i,4)));                           
        end  %gpu,����������������һЩ
    end
    
    APERTURE=ones(handles.green_Nrow, handles.green_Ncol);
    
    APER=[];
    handles.probe_ingreenNrow=[]; handles.probe_ingreenNcol=[]; handles.probestep=[];
    tic
    cuda_stem_potential;
    toc

    potential_temp = zeros(handles.green_Nrow,handles.green_Ncol, max(length(series_n), length(series_n_i)) );
    potential_temp(:) = potentialx(:) + sqrt(-1)*potentialy(:);
    
    for i=1:length(potential_temp(1,1,:)); %��Ҫת��
        potential(:,:,i) = potential_temp(:,:,i).';
    end
    clear potential_temp
end
%CPU�����Ǿ���
    %%��ɫ�ĳߴ磬��ɫ����ʸ��gx���ڼ���ԭ��λ�õģ�s2ʹ��peng�Ĳ�������ԭ���Ƴ��ֲ�


if get(handles.CPURB,'value')  %ʹ��CPU������ 
    
%����ÿ����Ƴ�
% potential=GetPotential4AllSlice_multicore_lobato_peng(handles.green_Ncol, handles.green_Nrow,... 
%     ele_n, absorp_n, ....�������ӻ���ԭ�ӵĵ��Ժ�����
%     ele_n_i, absorp_n_i, ... %ԭ�ӻ�ԭ��+���ӣ����Ի�������
%     series_n, series_n_i, ...  %ԭ�����д���   
%     s2_green, gx_green, gy_green, ...
%     sigma, PARAMETER, APERTURE, paraflag);   %Ϊ��STEM���㲻��������ֻ���뵽HRTEM��CBED��
APERTURE=ones(handles.green_Ncol, handles.green_Nrow);
potential=GetPotential4AllSlice_multicore_lobato_peng_corr(handles.green_Ncol, handles.green_Nrow,... 
    ele_n, absorp_n, ....�������ӻ���ԭ�ӵĵ��Ժ�����  %��������������peng�Ĳ���
    ele_n_i, absorp_n_i, ... %ԭ�ӻ�ԭ��+���ӣ����Ի�������
    series_n, series_n_i, ...  %ԭ�����д���   
    ele_n_corr, ele_n_i_corr, corr_info, ... %�������peng��ϵ��������������Ƴ�����
    series_n_corr, series_n_i_corr,...
    s2_green, gx_green, gy_green, ...
    sigma, PARAMETER, APERTURE, paraflag);   %Ϊ��STEM���㲻��������ֻ���뵽HRTEM��CBED��
end

%����ÿ����Ƴ� %ʹ�ö��̼߳���
% potential=GetPotential4AllSlice_multicore_lobato_peng(handles.green_Ncol, handles.green_Nrow,... 
%     ele_n, absorp_n, ....�������ӻ���ԭ�ӵĵ��Ժ�����
%     ele_n_i, absorp_n_i, ... %ԭ�ӻ�ԭ��+���ӣ����Ի�������
%     series_n, series_n_i, ...  %ԭ�����д���
%     s2_green, gx_green, gy_green, ...
%     sigma, PARAMETER, APERTURE, paraflag);   %Ϊ��STEM���㲻��������ֻ���뵽HRTEM��CBED��

% %����ÿ����Ƴ�
% potential=GetPotential4AllSlice(handles.green_Ncol, handles.green_Nrow,... 
%     ele_n, absorp_n, ....�������ӻ���ԭ�ӵĵ��Ժ�����
%     ele_n_i, absorp_n_i, ... %ԭ�ӻ�ԭ��+���ӣ����Ի�������
%     series_n, series_n_i, ...  %ԭ�����д���
%     s2_green, gx_green, gy_green, ...
%     sigma, PARAMETER, APERTURE);   %Ϊ��STEM���㲻��������ֻ���뵽HRTEM��CBED��

%�������
myap=ones(handles.CBEDprobesx);
myap(find(p2_probe>1/(9*handles.sampling*handles.sampling)))=0;  %16�Ǽ������� %9�Ǽ���2/3�Ĺ���

%����
myresul=zeros(handles.CBEDprobesx);
[probesx, probesy]=size(handles.probe);
for nn=1:length(handles.gfsf);
    mywave=handles.probe(:,:,nn);  %���벨����
    for k=1:length(potential(1,1,:));
        mywave =  mywave.*potential(handles.CBED_ingreenNrow:handles.CBED_ingreenNrow+handles.CBEDprobesx-1, ...
                   handles.CBED_ingreenNcol:handles.CBED_ingreenNcol+handles.CBEDprobesx-1,k);
        mywave=fft2(mywave);   %ע�⣬��ʡ��һ��fftshift
        mywave=mywave.*psf_fft;
        mywave=ifft2(mywave);
    end
    mywave=fftshift(fft2(mywave));
    myresul=myresul+handles.gfsf(nn)*(abs(mywave).^2);
end
%myresul=myresul.*myap;  %���Ϲ�����
return;

 

        
function [myresul,mid_ceng_mat] = STEMsimulation(handles, ele_n, ele_n_i, absorp_n, absorp_n_i, series_n, series_n_i, ele_n_corr, ele_n_i_corr, corr_info, ... %��������peng��ϵ��������������Ƴ�����
    series_n_corr, series_n_i_corr);
%��Ҫ����green�ĳߴ磬�Լ�����ԭ�ӵ����λ�ã�����������Լ�ƴ�ջ�ԭͼʱ�������
h_ba=(6.6255916e-34)/2/pi;
h=6.6255916e-34;
e=1.602102e-19;
me=9.1093807e-31;
c=2.9979251e+8;
%#############�����￪ʼ��gpu��cpu��ֻ�Է��ε���״���д�����������
    if handles.green_Nrow<handles.green_Ncol handles.green_Nrow=handles.green_Ncol;else handles.green_Ncol=handles.green_Nrow; end   %20201107
%#################################################################

PARAMETER=2*pi*(h_ba)*(h_ba)/(e*me*1.0E-10*1.0E-10)/(handles.green_Nrow*handles.green_Ncol*handles.sampling*handles.sampling);  %Angstrom��λ
sigma=2*pi/(handles.lambda*10*handles.vol/1000)*(me*c*c+e*handles.vol)/(2*me*c*c+e*handles.vol);  %��λ����һ����/kV*A
% %���ԣ���kirkland��һ��ͼ��ͼ5-2
% vol=10:10:1000;
% lambda=1.0e+9*h.*c./sqrt(e.*vol*1000.*(2*me*c*c+e.*vol*1000));
% sigma=2.*pi./(lambda.*10.*vol).*(me*c*c+e.*vol*1000)./(2*me*c*c+e.*vol*1000);

%APERTURE  % �������󣬱�֤�����������wrapЧӦ
APERTURE=ones(handles.green_Ncol, handles.green_Nrow);
[tx, ty]= meshgrid(-handles.green_Nrow/2:handles.green_Nrow/2-1, -handles.green_Ncol/2:handles.green_Ncol/2-1);
minn=min(handles.green_Nrow/2,handles.green_Ncol/2);
%APERTURE( find((tx./(handles.green_Nrow)).^2+(ty/(handles.green_Ncol)).^2>0.25/4) )=0;
clear tx ty

%%��ɫ�ĳߴ磬��ɫ����ʸ��gx���ڼ���ԭ��λ�õģ�s2ʹ��peng�Ĳ�������ԭ���Ƴ��ֲ�
[gx_green,gy_green]=meshgrid((-handles.green_Nrow/2:(handles.green_Nrow/2-1))./(handles.green_Nrow*handles.sampling), ...
    (-handles.green_Ncol/2:(handles.green_Ncol/2-1))./(handles.green_Ncol*handles.sampling)); %REPRO��λ��1/A)^2,����(1/nm)^2
sx_green=gx_green/2;  % peng ��s����
sy_green=gy_green/2;
s2_green=sx_green.^2+sy_green.^2;

% 
%   probe_propogation;��zheight�й� old codes20210119
% [px, py]=meshgrid( (-handles.probesx/2:(handles.probesx/2-1))./(handles.probesx*handles.sampling), ...
%     (-handles.probesx/2:(handles.probesx/2-1))./(handles.probesx*handles.sampling));
% p2_probe=(px.^2+py.^2);  %������
% psf_fft=exp(-sqrt(-1)*pi*handles.lambda*10*handles.eachthick*p2_probe);   %���ڵ�λA���㡣
% psf_fft(find(p2_probe>1/(16*handles.sampling*handles.sampling)))=0;  %��������
% psf_fft=ifftshift(psf_fft);
%add 20210119
Vol=str2num(get(handles.edit_Vol,'String')); %��ѹ
lambda=1.0e+9*h.*c./sqrt(e*Vol*1000*(2*me*c*c+e*Vol*1000));  %���㲨������λnm
tilt=str2num(get(handles.edit_Tilt_1,'String'))*10^(-3);
phitilt=-str2num(get(handles.edit_Tilt_2,'String'));  % ���ȵ�λ ���� 10mrad  
tiltx=tilt/lambda*cos(phitilt/180*pi); %���㵽���׸�ʽ�ĵ�λ 1/nm
tilty=tilt/lambda*sin(phitilt/180*pi);

%20210119 ������б�Ĳ������ڶ�㷨����ʱ��ע���----------------
%   ������propogation;��zheight�й�
[px, py]=meshgrid( (-handles.probesx/2:(handles.probesx/2-1))./(handles.probesx*handles.sampling), ...
    (-handles.probesx/2:(handles.probesx/2-1))./(handles.probesx*handles.sampling));
p2_probe=((px+tiltx*0.1).^2+(py+tilty*0.1).^2);  %������
psf_fft=exp(-sqrt(-1)*pi*handles.lambda*10*handles.eachthick*p2_probe);   %���ڵ�λA���㡣
%psf_fft(find(p2_probe>1/(16*handles.sampling*handles.sampling)))=0;  %��������
psf_fft(find(p2_probe>1/(9*handles.sampling*handles.sampling)))=0;  %����2/3�Ĺ���
psf_fft=ifftshift(psf_fft);

%end 20210119



aper_dect = str2num(get(handles.edit79,'string'));%������չ���
aper_range = aper_dect*0.001/(handles.lambda*10); %1/A ��λ
%�����м������չ���

%Kepeng
if get(handles.radiobutton11,'value')  %stem��
    num=length(aper_dect(:,1));
    APER=zeros(handles.probesx, handles.probesx,num);
    for i=1:num  %�������չ���
        temp=zeros(handles.probesx, handles.probesx);
        temp(find((px.^2+py.^2)>=min(aper_range(i,:)).^2 & (px.^2+py.^2)<=max(aper_range(i,:)).^2 )) = 1;
        APER(:,:,i)=temp;
    end
end

if get(handles.radiobutton13,'value')
    APER=makeAPER(handles);
end

% if 1  %3d stem
%     [dx,dy] = meshgrid( [1:handles.probesx]-mycen(handles.probesx) ,[1:handles.probesx]-mycen(handles.probesx) );
%     APER(:,:,end+1) = dx.*sum(APER,3); APER(:,:,end+1) = dy.*sum(APER(:,:,1:end-1),3);
% end
%Kepeng

% 
% myresul1 = cuda_STEM_core_complex(PARAMETER, sigma, ele_n, absorp_n, ele_n_i, absorp_n_i, series_n, series_n_i, gx_green,gy_green, s2_green,APERTURE, ... ...
%        handles.probe, handles.gfsf, psf_fft, handles.green_Ncol, handles.green_Nrow, handles.probe_ingreenNrow, handles.probe_ingreenNcol, ...
%        handles.width_red, handles.hight_red, handles.probestep,APER); 

paraflag = getparaflag(handles); %which parameter will be used for the scatting factor
if get(handles.CPURB, 'value') 
%     [myresul,mid_ceng_mat] = cuda_STEM_CPU_Mcore_lobato_peng(PARAMETER, sigma, ele_n, absorp_n, ele_n_i, absorp_n_i, series_n, series_n_i, gx_green,gy_green, s2_green,APERTURE, ... ...
%         handles.probe, handles.gfsf, psf_fft, handles.green_Ncol, handles.green_Nrow, handles.probe_ingreenNrow, handles.probe_ingreenNcol, ...
%         handles.width_red, handles.hight_red, handles.probestep,APER, paraflag, handles.mid_slice_num);
    [myresul,mid_ceng_mat] = cuda_STEM_CPU_Mcore_lobato_peng_corr(PARAMETER, sigma, ele_n, absorp_n, ele_n_i, absorp_n_i, series_n, series_n_i, ...
        ele_n_corr, ele_n_i_corr, corr_info, ... %��������peng��ϵ��������������Ƴ�����
        series_n_corr, series_n_i_corr,...
        gx_green,gy_green, s2_green,APERTURE, ... ...
        handles.probe, handles.gfsf, psf_fft, handles.green_Ncol, handles.green_Nrow, handles.probe_ingreenNrow, handles.probe_ingreenNcol, ...
        handles.width_red, handles.hight_red, handles.probestep,APER, paraflag, handles.mid_slice_num);
    
    
    % myresul1 = cuda_STEM_core_complex(PARAMETER, sigma, ele_n, absorp_n, ele_n_i, absorp_n_i, series_n, series_n_i, gx_green,gy_green, s2_green,APERTURE, ... ...
%        handles.probe, handles.gfsf, psf_fft, handles.green_Ncol, handles.green_Nrow, handles.probe_ingreenNrow, handles.probe_ingreenNcol, ...
%        handles.width_red, handles.hight_red, handles.probestep,APER); 
end
if get(handles.GPURB, 'value')  %���ʹ��GPU����
    if paraflag=='n'
        disp('Calculating the correction for Peng''s scattering factor')
        corr_info_matrix=0.*s2_green;%����ֵ
%         for i = 1:length(corr_info(:,1))  %��ÿ�ֵ������Ƴ����������������٣���������Ԫ�ظ����൱�ľ���ͺ�
%             r2 = s2_green;
%             r2(find(r2>=corr_info(i,5).^2 & r2<corr_info(i,6).^2))=(sin(pi * ((sqrt(r2(find(r2>=corr_info(i,5).^2 & r2<corr_info(i,6).^2)) )-corr_info(i,5))./(corr_info(i,6)-corr_info(i,5))-0.5) )+1)/2;
%             r2(find(r2<corr_info(i,5).^2))=0;
%             r2(find(r2>=corr_info(i,6).^2))=1;
%             corr_info_matrix(:,:,i) = r2.*( corr_info(i,1).*exp(-s2_green.*corr_info(i,3)) + corr_info(i,2).*exp(-s2_green.*corr_info(i,4)));                           
%         end  %gpu,����������������һЩ
        for i = 1:length(corr_info(:,1))  %��ÿ�ֵ������Ƴ����������������٣���������Ԫ�ظ����൱�ľ���ͺ�
            s2 = s2_green;
            g2 = 4*s2;
            corr_info_matrix(:,:,i) =(corr_info(i,11)*(2+corr_info(i,12)*g2)./(1+corr_info(i,12).*g2).^2 + ...
                                  corr_info(i,13)*(2+corr_info(i,14)*g2)./(1+corr_info(i,14).*g2).^2 + ...
                                  corr_info(i,15)*(2+corr_info(i,16)*g2)./(1+corr_info(i,16).*g2).^2 + ...
                                  corr_info(i,17)*(2+corr_info(i,18)*g2)./(1+corr_info(i,18).*g2).^2 + ...
                                  corr_info(i,19)*(2+corr_info(i,20)*g2)./(1+corr_info(i,20).*g2).^2 ) -...
                                 (corr_info(i,1).*exp(-s2.*corr_info(i,2)) ...
                                 + corr_info(i,3).*exp(-s2.*corr_info(i,4)) ...
                                 + corr_info(i,5).*exp(-s2.*corr_info(i,6)) ...
                                 + corr_info(i,7).*exp(-s2.*corr_info(i,8)) ...
                                 + corr_info(i,9).*exp(-s2.*corr_info(i,10)));                           
           end  %gpu,����������������һЩ
    end
    %cuda���㣬��Ҫ����ı����� corr_info_matrix; ele_n_corr, ele_n_i_corr;series_n_corr, series_n_i_corr
    tic
    cuda_stem;
    toc
end
 %save myresul myresul
 %stop;




% --- Executes on button press in pushbutton4.
function pushbutton4_Callback(hObject, eventdata, handles)
if get(handles.radiobutton11, 'value') | get(handles.radiobutton13, 'value')   %�����STEM��IDPC�ļ��㣬��ʾ��Ӧ������
    handles.sampling=str2num(get(handles.edit75,'string'));  %��/pixel   %����Ǵ��������ͶӰ�ƺ�������ĳ����ʣ�����probe��ɨ�������
    if ~isempty ( strfind(get(handles.edit77,'string'), '.'))
        msgbox('''Scan Step'' Must be an integer');
        return;
    end
    handles.probestepsampling=str2num(get(handles.edit77,'string')) * handles.sampling;  %��ʾɨ��ʱ��ķֱ���
    handles.probesx=str2num(get(handles.edit4,'string'));  %��ȡprobe�ĳߴ�
    if rem(handles.probesx,2)==1 handles.probesx=handles.probesx+1;end
    
    %�趨ɨ������topleft-rightdown������
        %��ɫ��Ŀ�Ⱥ͸߶ȷֱ�Ϊ��������
    width_red=str2num(get(handles.edit20,'string')); hight_red=str2num(get(handles.edit19,'string'));
    tl_red=[str2num(get(handles.edit17,'string')), str2num(get(handles.edit18,'string'))]; 
    rd_red=tl_red+[width_red, hight_red]*handles.probestepsampling;   %��λ��anstrong
    %���ϽǺ����½ǵ����꣨x��y��ֵ����λΪanstrong���ʼ��еĺ�ɫ��

    %���ڿ�ʼ��ʱ�򸳳�ֵoutside_ext=4;  %����������4����������֤ͶӰ�Ƴ�������Ҷ�任ʱ���������ص�Ч������Ӱ��probeɨ������
                %�ʼ��е���ɫ�򣬵�λҲ�ǰ�
    tl_green=tl_red-handles.sampling*round(handles.outside_ext./handles.sampling)-handles.probesx/2*handles.sampling;  %�ڶ��������Ĳ�������Ϊ�˱�֤��ɨ��ʱ����ѡ���������Ͻ�����������������
    rd_green=rd_red+handles.sampling*round(handles.outside_ext./handles.sampling)+handles.probesx/2*handles.sampling;  %�ڶ����ӷ��Ĳ�������Ϊ�˱�֤��ɨ��ʱ����ѡ���������Ͻ�����������������

    greensize = round((rd_green - tl_green)./handles.sampling);  
    if rem(greensize(1),2)==1 greensize(1)=greensize(1)+1;end %��ɫ�����ͼ��ߴ磬��λ������
    if rem(greensize(2),2)==1 greensize(2)=greensize(2)+1;end
      %����probeλ����ɫ���ڵ���ʼ�����λ��;�����ĳߴ��С��width_red��hight_red
    probe_in_green=[round(handles.outside_ext./handles.sampling)+1, round(handles.outside_ext./handles.sampling)+1];  

    %��Ҫ��һ����ɫ����������������ʣ�probe�ĳߴ磻��ɫ��Ŀ�͸ߡ�~~~~~~~~~~~~~
    load allresul
    atompos= allresul.x_new(:,2:4);  %�����޸�һ��ԭ�ӵ����꣬������Ϣ��Ҫ�鿴handles.x�����ȡ�Ľ����

    handles.green_Nrow=greensize(1);
    handles.green_Ncol=greensize(2);  %�к���
    handles.probe_ingreenNrow=probe_in_green(1);  
    handles.probe_ingreenNcol=probe_in_green(2);   %����ɫ����ѡȡ������λ��
    handles.width_red=width_red; 
    handles.hight_red=hight_red;  %ɨ���������ظ���
    handles.probestep=str2num(get(handles.edit77,'string')); %ɨ��ʱ�Ĳ���(��λ�Ǳ�����
    handles.tl_green=tl_green;   %����ɫ������Ͻ����꣨����λ��������������Ҫ֮�������е�ԭ��λ�ã�����ת��
    handles.rd_green=tl_green+greensize*handles.sampling;%��ɫ������½����꣬��Ҫ������Щԭ���������Χ��ģ��Ͳ�����stem����ĺ������㡣

    maxz=max(atompos(:,3));
    Drawsupercell_figure(allresul.x_new, handles.dis_atomsize, allresul.x_new(:,4), 0, 0);
    hold on; %plot(tl_green(1), tl_green(2), 'ko');plot(rd_green(1), rd_green(2), 'ko')
         line([tl_green(1),rd_green(1)],[tl_green(2),tl_green(2)],[maxz,maxz],'color','g','linewidth',3);
         line([tl_green(1),rd_green(1)],[rd_green(2),rd_green(2)],[maxz,maxz],'color','g','linewidth',3);
         line([tl_green(1),tl_green(1)],[tl_green(2),rd_green(2)],[maxz,maxz],'color','g','linewidth',3);
         line([rd_green(1),rd_green(1)],[tl_green(2),rd_green(2)],[maxz,maxz],'color','g','linewidth',3);
    hold on; %plot(tl_red(1), tl_red(2), 'ko');plot(rd_red(1), rd_red(2), 'ko')
         line([tl_red(1),rd_red(1)],[tl_red(2),tl_red(2)],[maxz,maxz],'color','r','linewidth',3);
         line([tl_red(1),rd_red(1)],[rd_red(2),rd_red(2)],[maxz,maxz],'color','r','linewidth',3);
         line([tl_red(1),tl_red(1)],[tl_red(2),rd_red(2)],[maxz,maxz],'color','r','linewidth',3);
         line([rd_red(1),rd_red(1)],[tl_red(2),rd_red(2)],[maxz,maxz],'color','r','linewidth',3);
         %��һ��probe����״�����probe������λ�ں�ɫ��������Ͻ�����
    tl_topleftprobe=tl_green+(probe_in_green-1)*handles.sampling; rd_topleftprobe=tl_green+(probe_in_green-1)*handles.sampling+handles.probesx*handles.sampling;
    hold on; %plot(tl_topleftprobe(1), tl_topleftprobe(2), 'ko');plot(rd_topleftprobe(1), rd_topleftprobe(2), 'ko')
         line([tl_topleftprobe(1),rd_topleftprobe(1)],[tl_topleftprobe(2),tl_topleftprobe(2)],[maxz,maxz],'color','k','linewidth',3);
         line([tl_topleftprobe(1),rd_topleftprobe(1)],[rd_topleftprobe(2),rd_topleftprobe(2)],[maxz,maxz],'color','k','linewidth',3);
         line([tl_topleftprobe(1),tl_topleftprobe(1)],[tl_topleftprobe(2),rd_topleftprobe(2)],[maxz,maxz],'color','k','linewidth',3);
         line([rd_topleftprobe(1),rd_topleftprobe(1)],[tl_topleftprobe(2),rd_topleftprobe(2)],[maxz,maxz],'color','k','linewidth',3);
    view([0,0,-1])
    axis equal
end

if get(handles.radiobutton10, 'value')  %�����CBEDͼ
    handles.sampling=str2num(get(handles.edit75,'string'));  %��/pixel   %����Ǵ��������ͶӰ�ƺ�������ĳ����ʣ�����probe��ɨ�������
    handles.CBEDprobesx=str2num(get(handles.edit4,'string'));  %��ȡprobe�ߴ�
    if rem(handles.CBEDprobesx,2)==1 handles.CBEDprobesx=handles.CBEDprobesx+1;end


    tl_red=[str2num(get(handles.edit17,'string')), str2num(get(handles.edit18,'string'))]; 
    %��ɫ������ͼ����������probe�ĳߴ�
    rd_red=tl_red+handles.CBEDprobesx*handles.sampling; 

    tl_green=tl_red-handles.sampling*round(handles.outside_ext./handles.sampling);  %����λ  %��ɫ�����ϽǺ����½����ꡣͼ����Ҫ��չһ���ߴ磬�������wrap
    rd_green=rd_red+handles.sampling*round(handles.outside_ext./handles.sampling);  

    greensize = round((rd_green - tl_green)./handles.sampling);  
    if rem(greensize(1),2)==1 greensize(1)=greensize(1)+1;end %��ɫ�����ͼ��ߴ磬��λ������
    if rem(greensize(2),2)==1 greensize(2)=greensize(2)+1;end
    %����imagingλ����ɫ���ڵ���ʼ�����λ��;֮��Ӵ�ͼ�вü����ճ�����
    imaging_in_green=[round(handles.outside_ext./handles.sampling)+1, round(handles.outside_ext./handles.sampling)+1];

    load allresul
    x=allresul.x_new;
    atompos=x(:,2:4);  % ��λ�䰣~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


    handles.atompos= atompos;  %�����޸�һ��ԭ�ӵ����꣬������Ϣ��Ҫ�鿴handles.x�����ȡ�Ľ����

    handles.green_Nrow=greensize(1);
    handles.green_Ncol=greensize(2);  %�к���
    handles.CBED_ingreenNrow=imaging_in_green(1); %ע���Ѿ�+1�ˣ����Դ���㿪ʼֱ��ȡλ�� 
    handles.CBED_ingreenNcol=imaging_in_green(2);   %����ɫ����ѡȡ������λ��
    handles.tl_green=tl_green;   %����ɫ������Ͻ����꣨����λ��������������Ҫ֮�������е�ԭ��λ�ã�����ת��
    handles.rd_green=tl_green+greensize*handles.sampling;%��ɫ������½����꣬��Ҫ������Щԭ���������Χ��ģ��Ͳ�����stem����ĺ������㡣

    maxz=max(atompos(:,3));
    load allresul
    Drawsupercell_figure(allresul.x_new, handles.dis_atomsize, allresul.x_new(:,4), 0, 0);
    hold on; %plot(tl_green(1), tl_green(2), 'ko');plot(rd_green(1), rd_green(2), 'ko')
         line([tl_green(1),rd_green(1)],[tl_green(2),tl_green(2)],[maxz,maxz],'color','g','linewidth',3);
         line([tl_green(1),rd_green(1)],[rd_green(2),rd_green(2)],[maxz,maxz],'color','g','linewidth',3);
         line([tl_green(1),tl_green(1)],[tl_green(2),rd_green(2)],[maxz,maxz],'color','g','linewidth',3);
         line([rd_green(1),rd_green(1)],[tl_green(2),rd_green(2)],[maxz,maxz],'color','g','linewidth',3);
    hold on; %plot(tl_red(1), tl_red(2), 'ko');plot(rd_red(1), rd_red(2), 'ko')
         line([tl_red(1),rd_red(1)],[tl_red(2),tl_red(2)],[maxz,maxz],'color','r','linewidth',3);
         line([tl_red(1),rd_red(1)],[rd_red(2),rd_red(2)],[maxz,maxz],'color','r','linewidth',3);
         line([tl_red(1),tl_red(1)],[tl_red(2),rd_red(2)],[maxz,maxz],'color','r','linewidth',3);
         line([rd_red(1),rd_red(1)],[tl_red(2),rd_red(2)],[maxz,maxz],'color','r','linewidth',3);
    view([0,0,-1])
    axis equal
end

if get(handles.radiobutton9, 'value')  %�����HRTEM�ļ��㣬��ʾ��Ӧ������
    handles.sampling=str2num(get(handles.edit75,'string'));  %��/pixel   %����Ǵ��������ͶӰ�ƺ�������ĳ����ʣ�����probe��ɨ�������
    handles.HRTEMedgesx=str2num(get(handles.edit4,'string'));  %��ȡ��Χ�ж��������Ҫ�õ�

    %�趨��������topleft-rightdown������
    %��ɫ��Ŀ�Ⱥ͸߶ȷֱ�Ϊ��������
    width_red=str2num(get(handles.edit20,'string')); hight_red=str2num(get(handles.edit19,'string'));
    
    tl_red=[str2num(get(handles.edit17,'string')), str2num(get(handles.edit18,'string'))]; 
    rd_red=tl_red+[width_red, hight_red]*handles.sampling;   %��λ��anstrong  %������STEM��ͬ
    %���ϽǺ����½ǵ����꣨x��y��ֵ����λΪanstrong��

    %��ɫ�򣬵�λҲ�ǰ�����ʵ�ʳ���ʱ�����ϽǺ����½����֮꣬��Ҫ�۵�����������������Ǳ߽���Եĳߴ磬����Ҫ����2��Ҳ����˵���۵�256�Ļ�������ʱ��ѡ��Ĵ�С����Ҫ512
    tl_green=tl_red-handles.HRTEMedgesx*handles.sampling;  %����λ
    rd_green=rd_red+handles.HRTEMedgesx*handles.sampling;  

    greensize = round((rd_green - tl_green)./handles.sampling);  
    if rem(greensize(1),2)==1 greensize(1)=greensize(1)+1;end %��ɫ�����ͼ��ߴ磬��λ������
    if rem(greensize(2),2)==1 greensize(2)=greensize(2)+1;end
    %����imagingλ����ɫ���ڵ���ʼ�����λ��;֮��Ӵ�ͼ�вü����ճ�����
    imaging_in_green=[handles.HRTEMedgesx+1, handles.HRTEMedgesx+1];  

    %��Ҫ��һ����ɫ����������������ʣ�probe�ĳߴ磻��ɫ��Ŀ�͸ߡ�~~~~~~~~~~~~~


    x=handles.x;
    atompos=x(:,2:4);  % ��λ�䰣~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    %��Ҫ����������е�ԭ��λ�ö�������ɫ�Ŀ��������ʼλ����ʲô

    handles.atompos= atompos;  %�����޸�һ��ԭ�ӵ����꣬������Ϣ��Ҫ�鿴handles.x�����ȡ�Ľ����

%     handles.eachthick=eachthick;  %��¼���
    handles.green_Nrow=greensize(1);
    handles.green_Ncol=greensize(2);  %�к���
    handles.HRTEM_ingreenNrow=imaging_in_green(1); %ע���Ѿ�+1�ˣ����Դ���㿪ʼֱ��ȡλ�� 
    handles.HRTEM_ingreenNcol=imaging_in_green(2);   %����ɫ����ѡȡ������λ��
    handles.tl_green=tl_green;   %����ɫ������Ͻ����꣨����λ��������������Ҫ֮�������е�ԭ��λ�ã�����ת��
    handles.rd_green=tl_green+greensize*handles.sampling;%��ɫ������½����꣬��Ҫ������Щԭ���������Χ��ģ��Ͳ�����stem����ĺ������㡣
    handles.width_red=width_red;
    handles.hight_red=hight_red;
    
    maxz=max(atompos(:,3));
    load allresul
    Drawsupercell_figure(allresul.x_new, handles.dis_atomsize, allresul.x_new(:,4), 0, 0);
    %axes(handles.axes1);
    hold on; %plot(tl_green(1), tl_green(2), 'ko');plot(rd_green(1), rd_green(2), 'ko')
         line([tl_green(1),rd_green(1)],[tl_green(2),tl_green(2)],[maxz,maxz],'color','g','linewidth',3);
         line([tl_green(1),rd_green(1)],[rd_green(2),rd_green(2)],[maxz,maxz],'color','g','linewidth',3);
         line([tl_green(1),tl_green(1)],[tl_green(2),rd_green(2)],[maxz,maxz],'color','g','linewidth',3);
         line([rd_green(1),rd_green(1)],[tl_green(2),rd_green(2)],[maxz,maxz],'color','g','linewidth',3);
    hold on; %plot(tl_red(1), tl_red(2), 'ko');plot(rd_red(1), rd_red(2), 'ko')
         line([tl_red(1),rd_red(1)],[tl_red(2),tl_red(2)],[maxz,maxz],'color','r','linewidth',3);
         line([tl_red(1),rd_red(1)],[rd_red(2),rd_red(2)],[maxz,maxz],'color','r','linewidth',3);
         line([tl_red(1),tl_red(1)],[tl_red(2),rd_red(2)],[maxz,maxz],'color','r','linewidth',3);
         line([rd_red(1),rd_red(1)],[tl_red(2),rd_red(2)],[maxz,maxz],'color','r','linewidth',3);
    view([0,0,-1])
    axis equal
end
guidata(hObject, handles);





function edit78_Callback(hObject, eventdata, handles)
function edit78_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function pushbutton4_CreateFcn(hObject, eventdata, handles)



function edit79_Callback(hObject, eventdata, handles)
%�󲨳�����mrad����ȷ����С��sampling rateӦ���Ƕ���
h_ba=(6.6255916e-34)/2/pi;
h=6.6255916e-34;
e=1.602102e-19;
me=9.1093807e-31;
c=2.9979251e+8;
vol=str2num(get(handles.edit_Vol,'string'));
lambda=1.0e+9*h.*c./sqrt(e*vol*1000*(2*me*c*c+e*vol*1000));  %���㲨������λnm

aper_dect = str2num(get(handles.edit79,'string'));  %detector �Ĺ����ߴ�
kmax=max(aper_dect(:))*10^(-3);

kk=kmax./(lambda*10);  %��ߵĵ��ռ��Ƶ���Ƕ��٣���λ��1/��

%�������ռ����С�����λΪ sampling��������ߵĵ��ռ�Ƶ��Ϊ1/sampling��
%���⣬����s�������ЧֵΪ6 ����֮һֵ������g��ֵ�Ļ���Ϊ2*s=12 ����֮һ��
%���Ǽ��㲨������ʱ�򣬻����1/2�Ĺ�������2/3�Ĺ������Ա�֤���ᷢ�� wrap �����
%���ԣ��Ƽ� 1/sampling=kk ҪС�� 6
%����6֮��������ǲ�׼�ģ����������ڷ���������Ϊ�˱�֤��Ƶ��Ϣ���У�����Ҫsamping rate�㹻С��

%������ߵĿռ�radֵΪ���١���Ϊ��6����Ч�ģ����Գ��Բ����󣬾�����ߵ�mradֵ
disp(strcat('Max aperture of detector is ', num2str(6*lambda*10 *1000),' mrad'));
disp('Because 6 1/A is maximum reciprocal-lattice-vector to calculate the projected potential')
%�����У�Ϊ�˷�ֹwrap��������̵��׿ռ䶼���˰��Ĺ�����
%Ϊ�˱�֤���˽ϸ߿ռ仹���е��ռ�ķ�����0.5/sampling ����ߵĵ��ռ�Ƶ�ʣ� 0.5*(0.5*1/samping)*lambda=alpha mrad
% %���� sampling=lambda*0.25/alpha ,ע�⻯������rad�ĵ�λ;����һ��0.5��-0.5/sampling ��0.5/sampling�����Ƶ��
% disp(strcat('Imaging sampling rate should be smaller than ', num2str(0.5*0.5*lambda*10/kmax),' A/pixel'));
%���� sampling=lambda*0.25/alpha ,ע�⻯������rad�ĵ�λ;����һ��0.5��-0.5/sampling
%��0.5/sampling�����Ƶ��;����ȡ����2/3�Ĺ���
disp(strcat('Imaging sampling rate should be smaller than ', num2str(0.5*2/3*lambda*10/kmax),' A/pixel'));
pp=1;



% --- Executes during object creation, after setting all properties.
function edit79_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in radiobutton10.
function radiobutton10_Callback(hObject, eventdata, handles)
if get(handles.radiobutton10,'value')==1
    set(handles.radiobutton9, 'value', 0);
    set(handles.radiobutton11, 'value', 0);
    set(handles.radiobutton13, 'value', 0);

end
if get(handles.radiobutton10,'value')==1
    set(handles.text115, 'visible', 'off');  %HRTEM��عر�
    set(handles.text106, 'visible', 'off');  %��STEM���
    set(handles.text108, 'visible','off');
    set(handles.text112, 'visible','off');
    set(handles.text110, 'visible','off');
    set(handles.edit77, 'visible','off');
    set(handles.edit79, 'visible','off');
    
    set(handles.text28, 'visible','off');
    set(handles.text103, 'visible','off');
    set(handles.edit20, 'visible','off');
    set(handles.edit19, 'visible','off');
    
    set(handles.text21, 'visible', 'off');   % aper��mrad��1/nm��λ
    set(handles.text125, 'visible', 'on');
    set(handles.text24, 'visible', 'on');   % convergence & aperture��mrad��1/nm��λ
    set(handles.text134, 'visible', 'off');
end
guidata(hObject, handles);


% --- Executes on button press in radiobutton9.
function radiobutton9_Callback(hObject, eventdata, handles)
if get(handles.radiobutton9,'value')==1
    set(handles.radiobutton10, 'value', 0);
    set(handles.radiobutton11, 'value', 0);
    set(handles.radiobutton13, 'value', 0);

end

if get(handles.radiobutton9,'value')==1
    set(handles.text115, 'visible', 'on');  %HRTEM���
    
    set(handles.text112, 'visible','on');
    set(handles.text110, 'visible','on');
    set(handles.edit79, 'visible','on');
    
    set(handles.text28, 'visible','on');
    set(handles.text103, 'visible','on');
    set(handles.edit20, 'visible','on');
    set(handles.edit19, 'visible','on');
    
   
    set(handles.text106, 'visible', 'off');  %��STEM���
    set(handles.text108, 'visible','off');
    set(handles.edit77, 'visible','off');
    
    set(handles.text112, 'visible','off');
    set(handles.text110, 'visible','off');
    set(handles.edit79, 'visible','off');
    
    set(handles.text21, 'visible', 'on');   % aper��mrad��1/nm��λ
    set(handles.text125, 'visible', 'off');
    set(handles.text134, 'visible', 'on');   % convergence & aperture��mrad��1/nm��λ
    set(handles.text24, 'visible', 'off');

else
    
    set(handles.text115, 'visible', 'off');  %HRTEM���

end

guidata(hObject, handles);



% --------------------------------------------------------------------
function Untitled_11_Callback(hObject, eventdata, handles)
%20200429�޸ģ���������ȫ����
probesx=str2num(get(handles.edit4,'String'));
   [para_part1, U, V]=CommonPara_TEM(handles, probesx, probesx);  %��CTF��PhasePlate�ù̶��Ĵ�С����
   para_part2 = readtccfromhandles_newTEM(hObject, handles, para_part1.lambda);   %ֻ��һ���������ѿ���е����ݶ���ȡ�ˣ�������֤�����ط���Ҫ��ȡ����ʱ����ͳһ�Ĵ���
    %������Ϊ�������֣�lambda��Ҫ����
   residual_aberr=WholeTCC_residual_newTEM(para_part2, U,V);  %����phase plate visualing the residual aberrations
    %-------------------------------------
   mytcc=WholeTCC2D_newTEM(para_part1, para_part2, U,V,2);
   
   myap=ones(probesx);  %�������������kmax�����涨��
   myap=myaperture(myap,1/(para_part1.sampling.*probesx),1/(para_part1.sampling.*probesx),0,para_part1.gmax*0.01/para_part1.lambda,0,0,1); %���10mrad���������1���ٳ���0.01��������rad���ٳ���nm��λ�Ĳ���
   num=length(mytcc(1,1,:));  %�����ж��ٸ�probe��Ҫɨ��ͼ��
   sum_probe=zeros(probesx);
   
   handles.gfsf=para_part1.gfsf;
   for i=1:num
        guiyihua=sum(sum(abs(ifft2(ifftshift(myap.*mytcc(:,:,i)))).^2));  %20201226 add one aper to normalize 
        guiyihua=sqrt(guiyihua);  %20201231��Ҫ�󿪸��ţ�������probe����3.59��3.68ʽ
        sum_probe = sum_probe + handles.gfsf(i)*myap.*mytcc(:,:,i)./guiyihua;
   end
   

  [sx,sy]=size(myap);
  temp = abs(fftshift(ifft2(ifftshift(sum_probe)))).^2;

figure;imshow(abs(fftshift(ifft2(ifftshift(sum_probe)))).^2,'XData',(-probesx/2:probesx/2-1).*para_part1.sampling,'YData',(-probesx/2:probesx/2-1).*para_part1.sampling,'DisplayRange',[]);axis on;xlabel('nm');
title('Intensity of Probe Wave')
figure;plot((-probesx/2:probesx/2-1).*para_part1.sampling, temp(sx/2+1,:) ,'r','linewidth', 2);
grid on
title('Intensity profile of Probe Wave')
figure;imshow(-angle(sum_probe.*myap),'XData',U(1, 1:end),'YData',V(1:end, 1),'DisplayRange',[]);axis on;xlabel('1/nm');
title('Phase of Probe Wave')


% --------------------------------------------------------------------
function Untitled_8_Callback(hObject, eventdata, handles)
CTFflag=0;
CTF_Callback(hObject, eventdata, handles, CTFflag)



function CTF_Callback(hObject, eventdata, handles, CTFflag)
sx=256;sy=256;
[para_part1, U, V]=CommonPara_TEM(handles, sx, sy);  %��CTF��PhasePlate�ù̶��Ĵ�С����
para_part2 = readtccfromhandles_newTEM(hObject, handles, para_part1.lambda);   %ֻ��һ���������ѿ���е����ݶ���ȡ�ˣ�������֤�����ط���Ҫ��ȡ����ʱ����ͳһ�Ĵ���
%������Ϊ�������֣�lambda��Ҫ����
residual_aberr=WholeTCC_residual_newTEM(para_part2, U,V);  %����phase plate visualing the residual aberrations
%-------------------------------------
%mytcc=WholeTCC2D(tccpara,defocus, U,V);
mytcc=WholeTCC2D_newTEM(para_part1, para_part2, U,V,1);
myresul=zeros(sx,sy);
for i=1:length(para_part1.gfsf);
     myresul=myresul+para_part1.gfsf(i).*mytcc(:,:,i);
end

figure;hold off;
if CTFflag==0;
    myresul=real(myresul);
elseif CTFflag==1
    myresul=imag(myresul);
elseif CTFflag==2
    myresul=abs(myresul);
end
myap=ones(size(U));
myap=myaperture(myap,1/(para_part1.sampling.*sx),1/(para_part1.sampling.*sy),0,para_part1.gmax,0,0,1);
disimage=myresul.*myap;

imshow(disimage,'XData',U(1, 1:end).*para_part2.lambda*1000,'YData',V(1:end, 1).*para_part2.lambda*1000,'DisplayRange',[]);axis on;xlabel('mrad');

colorbar;
if CTFflag==0;
    title('Real CTF');
elseif CTFflag==1
    title('Imaginary CTF');
elseif CTFflag==2
    title('Damping CTF');
end
return;


% --- Executes on button press in radiobutton11.
function radiobutton11_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton11 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if get(handles.radiobutton11,'value')==1 
    set(handles.radiobutton10, 'value', 0);
    set(handles.radiobutton9, 'value', 0);
    set(handles.radiobutton13, 'value', 0);
end
if get(handles.radiobutton11,'value')==1
    set(handles.text115, 'visible', 'off');  %HRTEM��عر�
    
    set(handles.text21, 'visible', 'off');   % aper��mrad��1/nm��λ
    set(handles.text125, 'visible', 'on');
    set(handles.text24, 'visible', 'on');   % convergence & aperture��mrad��1/nm��λ
    set(handles.text134, 'visible', 'off');
    
    set(handles.text106, 'visible', 'on');  %STEM���
    set(handles.text108, 'visible','on');
    set(handles.text112, 'visible','on');
    set(handles.text110, 'visible','on');
    set(handles.edit77, 'visible','on');
    set(handles.edit79, 'visible','on');
    
    set(handles.text28, 'visible','on');
    set(handles.text103, 'visible','on');
    set(handles.edit20, 'visible','on');
    set(handles.edit19, 'visible','on');
else
    set(handles.text106, 'visible', 'off');
    set(handles.text108, 'visible','off');
    set(handles.text112, 'visible','off');
    set(handles.text110, 'visible','off');
    set(handles.edit77, 'visible','off');
    set(handles.edit79, 'visible','off');
    
    set(handles.text28, 'visible','off');
    set(handles.text103, 'visible','off');
    set(handles.edit20, 'visible','off');
    set(handles.edit19, 'visible','off');
    
    set(handles.text21, 'visible', 'on');   % aper��mrad��1/nm��λ
    set(handles.text125, 'visible', 'off');
    set(handles.text24, 'visible', 'off');   % convergence & aperture��mrad��1/nm��λ
    set(handles.text134, 'visible', 'on');
    
end
guidata(hObject, handles);
% Hint: get(hObject,'Value') returns toggle state of radiobutton11


% --------------------------------------------------------------------
function Untitled_13_Callback(hObject, eventdata, handles)
sx=256; sy=256;
[para_part1, U, V]=CommonPara_TEM(handles, sx, sy);  %��CTF��PhasePlate�ù̶��Ĵ�С����
para_part2 = readtccfromhandles_newTEM(hObject, handles, para_part1.lambda);   %ֻ��һ���������ѿ���е����ݶ���ȡ�ˣ�������֤�����ط���Ҫ��ȡ����ʱ����ͳһ�Ĵ���
%������Ϊ�������֣�lambda��Ҫ����

residual_aberr=WholeTCC_residual_newTEM(para_part2, U,V);  %����phase plate visualing the residual aberrations
residual_phase=angle(exp(sqrt(-1)*2*pi*residual_aberr/para_part1.lambda));


figure; hold off; imshow(residual_phase,'XData',U(1, 1:end).*para_part2.lambda*1000,'YData',V(1:end, 1).*para_part2.lambda*1000,'DisplayRange',[]);axis on;xlabel('mrad');title('ctg(kai)');
axis equal
    
prompt = {'Radiu 1'; 'Radiu 2'};
dlg_title = '2D Phase Plate';
num_lines = 1;
def = {'16';'28'};
answer = inputdlg(prompt,dlg_title,num_lines,def);
if isempty(answer) 
    return; 
end
rad1=str2num(answer{1}); rad2=str2num(answer{2}); %�����������Ѳ�����ɢ��ͼ��Ȧ����

title(strcat('2D Phase Plate ',num2str(rad1), '&',num2str(rad2)));
colorbar
hold on; grid on  %�����mrad��Բ����
x=linspace(-rad1,rad1,sx*5);
y=sqrt(rad1.*rad1-x.*x);
plot(x,y,'r-','LineWidth',3)
plot(x,-y,'r-','LineWidth',3)
%rad2=0.5*rad1;  %��1/2��aper����ֵ����ͼ��ֻ�ǻ��㵽mrad��λ
x=linspace(-rad2,rad2,sx*20);
y=sqrt(rad2.*rad2-x.*x);
plot(x,y,'b-','LineWidth',3)
plot(x,-y,'b-','LineWidth',3)
%title(strcat(num2str(rad1),' and  ', num2str(rad2), 'mrad phase plate due to the residual aberrations'))

guidata(hObject, handles);


% --------------------------------------------------------------------
function Untitled_12_Callback(hObject, eventdata, handles)
CTFflag=1;
CTF_Callback(hObject, eventdata, handles, CTFflag);
pp=1;


% --------------------------------------------------------------------
function Untitled_14_Callback(hObject, eventdata, handles)
CTFflag=2;
CTF_Callback(hObject, eventdata, handles, CTFflag);



function [tccpara, U, V]=CommonPara_TEM(handles, sx, sy);
h_ba=(6.6255916e-34)/2/pi;
h=6.6255916e-34;
e=1.602102e-19;
me=9.1093807e-31;
c=2.9979251e+8;
tccpara.sampling=str2num(get(handles.edit75,'String'))*0.1; %  ������ %UNIT NM/PIXEL  %������EW_RECONSTRUCT��һ���ĵ�λ
tccpara.mtotal=(get(handles.popupmenu4,'Value')-1)*2+1; %  �ܹ��ж��ٸ��뽹ƫ�룬�������ֵΪ5���ͱ�����1 2 3 4 5��mytcc�������е�����ͼ���뽹ƫ��

tccpara.alafa=str2num(get(handles.edit_Con,'String'))*(10^(-3)); % ��ǿ��
Vol=str2num(get(handles.edit_Vol,'String')); %��ѹ
tccpara.lambda=1.0e+9*h.*c./sqrt(e*Vol*1000*(2*me*c*c+e*Vol*1000));  %���㲨������λnm
strcat( 'wavelength of incident electron (nm): ', num2str(tccpara.lambda) )

tcccpara.yita=str2num(get(handles.edit_Spr,'String')); 
tccpara.gmax=str2num(get(handles.edit_Ape,'String'))*0.1;   %���㵽A�ĵ�λ
 
resul=Gaussian_focal(tccpara.mtotal,tcccpara.yita,tccpara.lambda,tccpara.gmax); %Ϊ�뽹ͼ����ռ��ͼ��İٷ�֮��
tccpara.gfsf=resul.gfsf;
tccpara.delta_yita=resul.delta_yita;

%����б���ٽǶ�
tccpara.tilt=str2num(get(handles.edit_Tilt_1,'String'))*10^(-3);tccpara.phitilt=-str2num(get(handles.edit_Tilt_2,'String'));  % ���ȵ�λ ���� 10mrad  
tccpara.tiltx=tccpara.tilt/tccpara.lambda*cos(tccpara.phitilt/180*pi); %���㵽���׸�ʽ�ĵ�λ 1/nm
tccpara.tilty=tccpara.tilt/tccpara.lambda*sin(tccpara.phitilt/180*pi);


v=-round((sx+1)/2)+1:sx-(round((sx+1)/2));
u=-round((sy+1)/2)+1:sy-(round((sy+1)/2));
u=u.*1/(tccpara.sampling.*sy);
v=v.*1/(tccpara.sampling.*sx);
[U,V]=meshgrid(u,v);
U=U+tccpara.tiltx;
V=V+tccpara.tilty;


function tccpara=readtccfromhandles_newTEM(hObject, handles, lambda)  %�����Ǹ���TEM�Ĳ�������ȡ��
tccpara.lambda=lambda;
%if get(handles.checkbox_polar,'Value')==1  %�ض��Ǽ�����ϵ
   %��ȡ˵����%�ο�Phil. Trans. R. Soc. A 2009, vol
%367:3755-3771;����Ȧ������ϵ����UM  1998 72 PP109-119�Լ�UM 1996 64 249-264�У�B2 S3 D4 B4
%S5 R5(D5)����û��ϵ������ˣ���ȡ���ݺ�Ӧ��Ҫ�ֱ���� 3  4  5  5  6  6

   tccpara.focus=str2num(get(handles.edit_Focus,'String'));
   
   %A1����ֵ���䣻A1�ĽǶ��ǵ羵��ʾ����ֵ��1/2 
   tccpara.A1=str2num(get(handles.edit_A1_1,'String')); tccpara.phiA1=-str2num(get(handles.edit_A1_2,'String'));
   
   %A2����ֵ���䣻A2�ĽǶ��ǵ羵��ʾ����ֵ��1/3
   tccpara.A2=str2num(get(handles.edit_A2_1,'String')); tccpara.phiA2=-str2num(get(handles.edit_A2_2,'String'));
  
   %B2����ֵ���䣻B2�ĽǶ���羵��ʾ��ȡ����
   tccpara.B2=str2num(get(handles.edit_B2_1,'String')); %�����������ʱ���rew�����3��������rew����Ϊ2000�����������6000����Ϊ�����е�ϵ�������1/3
   tccpara.phiB2=-str2num(get(handles.edit_B2_2,'String'));  
  
   %Cs�ĵ�λ��um��nm
   tccpara.Cs=str2num(get(handles.edit_Cs,'String'))*(10^3); %UNIT from um to nm
   
   %A3�ĵ�λ��um��nm��A3�ĽǶ���羵��ʾ�ĳ���4
   tccpara.A3=str2num(get(handles.edit_A3_1,'String'))*(10^3); %UNIT from um to nm
   tccpara.phiA3=-str2num(get(handles.edit_A3_2,'String'));
   
   %S3�ĵ�λ��um��nm��S3�ĽǶ���羵��ʾ��ȡ��������2
   tccpara.S3=str2num(get(handles.edit_S3_1,'String'))*(10^3); 
   tccpara.phiS3=-str2num(get(handles.edit_S3_2,'String'));
   
   %A4�ĵ�λ��um��nm��A4�ĽǶ���羵��ʾ�ĳ���5
   tccpara.A4=str2num(get(handles.edit_A4_1,'String'))*(10^3);
   tccpara.phiA4=-str2num(get(handles.edit_A4_2,'String'));
   
   %D4�ĵ�λ��um��nm��D4�ĽǶ���羵��ʾ��ȡ������3
   tccpara.D4=str2num(get(handles.edit_D4_1,'String'))*(10^3); 
   tccpara.phiD4=-str2num(get(handles.edit_D4_2,'String'));
   
   %B4�ĵ�λ��um��nm��B4�ĽǶ���羵��ʾ��һ��   %��Ҫ���о��ĽǶȹ�ϵ������
   tccpara.B4=str2num(get(handles.edit_B4_1,'String'))*(10^3); 
   tccpara.phiB4=-str2num(get(handles.edit_B4_2,'String'));
   
    %A5�ĵ�λ��mm��nm��A4�ĽǶ���羵��ʾ�ĳ���6
   tccpara.A5=str2num(get(handles.edit_A5_1,'String'))*(10^6); 
   tccpara.phiA5=-str2num(get(handles.edit_A5_2,'String'));
   
   %A5�ĵ�λ��mm��nm��
   tccpara.C5=str2num(get(handles.edit_C5,'String'))*(10^6);
   
   %S5�ĵ�λ��mm��nm�� %�����ֵ��hrtem���滹û�г��֣���˽Ƕȹ�ϵ��ʱû�й�
   tccpara.S5=str2num(get(handles.edit_S5_1,'String'))*(10^6); tccpara.phiS5=-str2num(get(handles.edit_S5_2,'String'));
  
   %D5�ĵ�λ��mm��nm��%�����ֵ��hrtem���滹û�г��֣���˽Ƕȹ�ϵ��ʱû�й�
   tccpara.D5=str2num(get(handles.edit_D5_1,'String'))*(10^6); tccpara.phiD5=-str2num(get(handles.edit_D5_2,'String'));

   
   
   tccpara.A1x=tccpara.A1*cos(tccpara.phiA1/180*pi);   tccpara.A1y=tccpara.A1*sin(tccpara.phiA1/180*pi);
   %C1=focus;  focus���������defocs������ƫ���
   tccpara.A2x=tccpara.A2*cos(tccpara.phiA2/180*pi);   tccpara.A2y=tccpara.A2*sin(tccpara.phiA2/180*pi);
   tccpara.B2x=tccpara.B2*cos(tccpara.phiB2/180*pi);     tccpara.B2y=tccpara.B2*sin(tccpara.phiB2/180*pi);
   tccpara.A3x=tccpara.A3*cos(tccpara.phiA3/180*pi);   tccpara.A3y=tccpara.A3*sin(tccpara.phiA3/180*pi);
   tccpara.S3x=tccpara.S3*cos(tccpara.phiS3/180*pi);   tccpara.S3y=tccpara.S3*sin(tccpara.phiS3/180*pi);
   tccpara.C3=tccpara.Cs;
   tccpara.A4x=tccpara.A4*cos(tccpara.phiA4/180*pi);   tccpara.A4y=tccpara.A4*sin(tccpara.phiA4/180*pi);
   tccpara.D4x=tccpara.D4*cos(tccpara.phiD4/180*pi);   tccpara.D4y=tccpara.D4*sin(tccpara.phiD4/180*pi);
   tccpara.B4x=tccpara.B4*cos(tccpara.phiB4/180*pi);     tccpara.B4y=tccpara.B4*sin(tccpara.phiB4/180*pi);
   tccpara.A5x=tccpara.A5*cos(tccpara.phiA5/180*pi);   tccpara.A5y=tccpara.A5*sin(tccpara.phiA5/180*pi);
   tccpara.S5x=tccpara.S5*cos(tccpara.phiS5/180*pi);   tccpara.S5y=tccpara.S5*sin(tccpara.phiS5/180*pi);
   tccpara.C5=tccpara.C5;
   tccpara.D5x=tccpara.D5*cos(tccpara.phiD5/180*pi);   tccpara.D5y=tccpara.D5*sin(tccpara.phiD5/180*pi);

pp=1;

function W=WholeTCC_residual_newTEM(para,u,v);   %���չ�ʽ����TCC������ϵ�������ǽ�ȥ
%��� formular��ʽ���ù�ʽ�����еĸ߽���ɢ��ʽ�����ֳ�����������˶�u��v������ݶ�

oumigau=u.*para.lambda;  %
oumigav=v.*para.lambda;
i=sqrt(-1);
%'2-fold astigmatism A1'
A1x=para.A1x;   A1y=para.A1y;
%C1=focus;  focus���������defocs������ƫ���
A2x=para.A2x;   A2y=para.A2y;
B2x=para.B2x;   B2y=para.B2y;
A3x=para.A3x;   A3y=para.A3y;
S3x=para.S3x;   S3y=para.S3y;
C3=para.Cs;
A4x=para.A4x;   A4y=para.A4y;
D4x=para.D4x;   D4y=para.D4y;
B4x=para.B4x;   B4y=para.B4y;
A5x=para.A5x;   A5y=para.A5y;
S5x=para.S5x;   S5y=para.S5y;
C5=para.C5;
D5x=para.D5x;   D5y=para.D5y;

%�ο�Phil. Trans. R. Soc. A 2009, vol
%367:3755-3771;��ʽ2.4ʽ�ӣ�oumiga=u*lambda+i*v*lambda;  ���Ultramicroscopy 64
%1996, pp.249-2264��equation 11
W= real( ((A2x/3 + (A2y*i)/3).*(oumigau - i*oumigav).^3) ) ...
    +real (((B2x + (B2y*i)).*(oumigau + i*oumigav).^2.*(oumigau - i*oumigav)) ) ...
    +real ( ((A3x/4 + (A3y*i)/4).*(oumigau - i*oumigav).^4) ) ...
    +real ( ((oumigau - i*oumigav).*(oumigau + i*oumigav).^3.*(S3x + (S3y*i)))  ) ...
    +  real(C3*(oumigau + i*oumigav).^2.*(oumigau - i*oumigav).^2)/4 ...
    + real( ((A4x/5 + (A4y*i)/5).*(oumigau - i*oumigav).^5) ) ...
    + real( ((D4x + (D4y*i)).*(oumigau + i*oumigav).^4.*(oumigau - i*oumigav)) ) ...
    + real( ((B4x + (B4y*i)).*((oumigau + i*oumigav).^3).*((oumigau - i*oumigav).^2)))  ...
    + real( ((A5x/6 + (A5y*i)/6).*(oumigau - i*oumigav).^6) ) ...
    + real( (oumigau + i*oumigav).^4.*(oumigau - i*oumigav).^2.*(S5x + (S5y*i)) ) ...
    + real(C5*(oumigau + i*oumigav).^3.*(oumigau - i*oumigav).^3)/6 ...
    + real( ((D5x + (D5y*i)).*(oumigau + i*oumigav).^5.*(oumigau - i*oumigav)) )  ;

    
    pp=1;


function [tlittle_all, damping]=WholeTCC2D_newTEM(para1, para2, u,v, flag);   %���չ�ʽ����TCC������ϵ�������ǽ�ȥ
%��� formular��ʽ���ù�ʽ�����еĸ߽���ɢ��ʽ�����ֳ�����������˶�u��v������ݶ�

oumigau=u.*para1.lambda;  %
oumigav=v.*para1.lambda;
i=sqrt(-1);
%'2-fold astigmatism A1'
A1x=para2.A1x;   A1y=para2.A1y;
C1=para2.focus;
A2x=para2.A2x;   A2y=para2.A2y;
B2x=para2.B2x;   B2y=para2.B2y;
A3x=para2.A3x;   A3y=para2.A3y;
S3x=para2.S3x;   S3y=para2.S3y;
C3=para2.Cs;
A4x=para2.A4x;   A4y=para2.A4y;
D4x=para2.D4x;   D4y=para2.D4y;
B4x=para2.B4x;   B4y=para2.B4y;
A5x=para2.A5x;   A5y=para2.A5y;
S5x=para2.S5x;   S5y=para2.S5y;
C5=para2.C5;
D5x=para2.D5x;   D5y=para2.D5y;

%�ο�Phil. Trans. R. Soc. A 2009, vol
%367:3755-3771;��ʽ2.4ʽ�ӣ�oumiga=u*lambda+i*v*lambda;
kai=WholeTCC_residual_newTEM(para2,u,v);
W=real(((A1x + (A1y*i)).*(oumigau - i*oumigav).^2))/2 ...   %'2-fold astigmatism A1'
    + real(C1.*(oumigau + i*oumigav).*(oumigau - i*oumigav))/2 + kai;

flag=2;
if flag==2  %��ԭ������ȣ�B�ȶ��б仯
%new
tempgradient_kai_u=para1.lambda.*(A1x*u+A1y*v) ...%    +C1.*para.lambda.*u ...
    +A2x*para1.lambda.^2/3*(3*u.*u-3*v.*v)+A2y*para1.lambda.^2/3*6*u.*v ...
    +B2x.*(3*u.^2+v.^2).*para1.lambda.^2 - B2y.*para1.lambda.^2.*2.*u.*v ...
    +1/4*para1.lambda.^3 * (A3x*(4*u.^3-12*u.*v.*v)+A3y*(12*u.*u.*v-4*v.^3)) ...
    +para1.lambda.^3* (2*u.*(S3x*(u.*u-v.*v)-S3y*2*u.*v) +(u.*u+v.*v).*(2*S3x*u-2*S3y*v)) ...
    +C3*para1.lambda^3*(u.^2+v.^2).*u ...
    +1/5*para1.lambda^4 * (A4x*(5*u.^4+5*v.^4-30.*u.^2.*v.^2)+A4y*(20*u.^3.*v-20*u.*v.^3)) ...
    +para1.lambda^4 * ( (2*u) .* (D4x*(u.^3-3*u.*v.^2)-D4y*(3*u.^2.*v-v.^3)) + (u.^2+v.^2).*(D4x*(3.*u.^2-3*v.^2)-D4y*6*u.*v) ) ...
    +para1.lambda^4 * (  2*(u.^2+v.^2).*2.*u .*(B4x*u-B4y*v) + (u.^2+v.^2).^2*B4x) ...
    +1/6*para1.lambda^5 * (A5x*(4.*u.^3-14*v).*(u.^2-v.^2)+A5x*(u.^4+v.^4-14*u.*v)*2.*u   - A5y*2*v.*(3*u.^4+3*v.^4-10*u.^2.*v.^2)-A5y*2.*v.*u.*(12*u.^3-20*u.*v.^2) ) ...
    +para1.lambda^5 * (2*(u.^2+v.^2).*2.*u.*(S5x*(u.^2-v.^2)-2*u.*v*S5y) + (u.^2+v.^2).^2.*(2*S5x*u-2*v*S5y) ) ...
    +1/6*para1.lambda^5 * C5 * 3*(u.^2+v.^2).^2.*2.*u ...
    +para1.lambda^5 * (2*u.*(D5x*(u.^4+v.^4-6*u.^2.*v.^2)-D5y*4*u.*v.*(u.^2-v.^2)) + (u.^2+v.^2).*(D5x.*(4.*u.^3-12*u.*v.^2)-D5y*(12*u.^2.*v-4*v.^3)) ) ;
    
    
tempgradient_kai_v=para1.lambda.*(A1y*u-A1x*v) ...%    +C1.*para.lambda.*v ...
    -A3x.*para1.lambda.^2/3*6*u.*v + A3y*para1.lambda.^2/3*(-3*v.^2+3*u.^2) ...
    + B2x.*para1.lambda.*para1.lambda.*2*u.*v - B2y.*para1.lambda.*para1.lambda.*(u.^2+3*v.^2) ...
    +1/4*para1.lambda.^3 * (A3x*(4*v.^3-12*u.*u.*v)+A3y*(4*u.^3-12*u.*v.*v)) ...
    +para1.lambda.^3* (2*v.*(S3x*(u.*u-v.*v)-S3y*2*u.*v) +(u.*u+v.*v).*(-2*S3x*v-2*S3y*u) ) ...
    +C3*para1.lambda.^3.*(u.^2+v.^2).*v ...
    +1/5*para1.lambda^4 * (A4x*(20*u.*v.^3-20.*u.^3.*v)+A4y*(5*v.^4+5*u.^4-30*u.^2.*v.^2)) ...
    +para1.lambda^4 * ( (2*v) .* (D4x*(u.^3-3*u.*v.^2)-D4y*(3*u.^2.*v-v.^3)) + (u.^2+v.^2).*(-D4x*6*u.*v-D4y*(3*u.^2-3*v.^2))) ...
    +para1.lambda^4 * (  2*(u.^2+v.^2).*2.*v .*(B4x*u-B4y*v) - (u.^2+v.^2).^2*B4y) ...
    +1/6*para1.lambda^5 * (A5x*(4.*v.^3-14*u).*(u.^2-v.^2)+A5x*(u.^4+v.^4-14*u.*v).*(-2.*v)   - A5y*2*u.*(3*u.^4+3*v.^4-10*u.^2.*v.^2)-A5y*2.*v.*u.*(12*v.^3-20*u.^2.*v) ) ...
    +para1.lambda^5 * (2*(u.^2+v.^2).*2.*v.*(S5x*(u.^2-v.^2)-2*u.*v*S5y) + (u.^2+v.^2).^2.*(-2*S5x*v-2*u*S5y) ) ...
    +1/6*para1.lambda^5 * C5 * 3*(u.^2+v.^2).^2.*2.*v ...
    +para1.lambda^5 * (2*v.*(D5x*(u.^4+v.^4-6*u.^2.*v.^2)-D5y*4*u.*v.*(u.^2-v.^2)) + (u.^2+v.^2).*(D5x.*(4.*v.^3-12*u.^2.*v)-D5y*(4*u.^3-12*u.*v.^2)) ) ;

end
if flag==1
    tempgradient_kai_u=0;
    tempgradient_kai_v=0;
end

for m=1:para1.mtotal
%    C1=focus+(m-(para.mtotal+1)/2)*para.delta_yita;
    
    gradient_kai_u=tempgradient_kai_u+C1.*para1.lambda.*u;
    gradient_kai_v=tempgradient_kai_v+C1.*para1.lambda.*v;

%gradient_W=gradient_W_u + i* gradient_W_v;
%�ο�um 1996 vol 64:109-135����ʽ5
E_s_coh=exp(-para1.alafa.*para1.alafa.*pi*pi/para1.lambda/para1.lambda*(gradient_kai_u.^2+gradient_kai_v.^2));
kai=2*pi*W/para1.lambda;
%�ο�Phil. Trans. R. Soc. A 2009, vol 367:3755-3771;��ʽ2.2�·�
%E_s_coh=exp(-para.alpha.*para.alpha/4/para.lambda/para.lambda.*gradient_W.*conj(gradient_W).*(2*pi*2*pi/para.lambda/para.lambda));  %�ռ������ ����para.alphaΪ���;gradient_kaiʵ���������ϵ�gradient_W����Ҫ�ٳ���w*pi��
%����ʽ����Phil. Trans. R. Soc. A 2009, vol
%367:3755-3771;��ʽ2.1д����ע��gradient_W_u_vʵ�����Ƕ�W�󵼵ġ�
%��ο�um 1996 vol 64:109-135��ȣ��������е� kai=W/lambda��WΪphil�����ױ�����
%���Բο�um 1996 vol 64:109-135����ʽ5,�뱾��ʾ���һ��
E_f_coh=exp(-i*pi*(m-(para1.mtotal+1)/2)*para1.delta_yita*(u.^2+v.^2).*para1.lambda);  %ɫ��   �����u��v���ǵ���ʸ��
%�ο�um 1996 vol 64:109-135����ʽ35,�뱾��ʽ��Դ
t=exp(-i*kai);
% if m==1
% figure; imshow(real(t),[]);title('real t')
% end
% figure; imshow(real(E_f_coh),[]);title('E f coh real')
% figure; imshow(E_s_coh,[]);title('E s coh')
%�ο�Phil. Trans. R. Soc. A 2009, vol 367:3755-3771;��ʽ2.2�·�
tlittle_all(:,:,m)=E_s_coh.*E_f_coh.*t;
damping(:,:,m)=E_s_coh.*E_f_coh;
%�˲�
[sx,sy]=size(u);
if min(sx,sy)>2  %�������һά���ݣ�������ͨ�˲�
%   tlittle_all(:,:,m)=myaperture(tlittle_all(:,:,m),1/(para.sampling.*sx),1/(para.sampling.*sy),0,para.gmax,0,0,1);
else 
    
   tlittle_all(1,round(para1.gmax./para1.sampling):end,m)=0;  %һά�����Ĵ�ͨ����Ƶ����Ϊ0����ʱpara.sampling�Ǳ�ʾ���ռ�ļ��Ƶ����0.01nm_1��ֵ���̶�������
end
end

   
function W=WholeTCC2D_newTEMphase(para1, para2, u,v);   %���չ�ʽ����TCC������ϵ�������ǽ�ȥ
%��� formular��ʽ���ù�ʽ�����еĸ߽���ɢ��ʽ�����ֳ�����������˶�u��v������ݶ�

oumigau=u.*para1.lambda;  %
oumigav=v.*para1.lambda;
i=sqrt(-1);
%'2-fold astigmatism A1'
A1x=para2.A1x;   A1y=para2.A1y;
C1=para2.focus;
A2x=para2.A2x;   A2y=para2.A2y;
B2x=para2.B2x;   B2y=para2.B2y;
A3x=para2.A3x;   A3y=para2.A3y;
S3x=para2.S3x;   S3y=para2.S3y;
C3=para2.Cs;
A4x=para2.A4x;   A4y=para2.A4y;
D4x=para2.D4x;   D4y=para2.D4y;
B4x=para2.B4x;   B4y=para2.B4y;
A5x=para2.A5x;   A5y=para2.A5y;
S5x=para2.S5x;   S5y=para2.S5y;
C5=para2.C5;
D5x=para2.D5x;   D5y=para2.D5y;

%�ο�Phil. Trans. R. Soc. A 2009, vol
%367:3755-3771;��ʽ2.4ʽ�ӣ�oumiga=u*lambda+i*v*lambda;
kai=WholeTCC_residual_newTEM(para2,u,v);
W=real(((A1x + (A1y*i)).*(oumigau - i*oumigav).^2))/2 ...   %'2-fold astigmatism A1'
    + real(C1.*(oumigau + i*oumigav).*(oumigau - i*oumigav))/2 + kai;

% 
% function edit83_Callback(hObject, eventdata, handles)
% function edit83_CreateFcn(hObject, eventdata, handles)
% if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
%     set(hObject,'BackgroundColor','white');
% end

% 
% 
% function edit84_Callback(hObject, eventdata, handles)
% function edit84_CreateFcn(hObject, eventdata, handles)
% if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
%     set(hObject,'BackgroundColor','white');
% end






function edit92_Callback(hObject, eventdata, handles)
function edit92_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end




% --- Executes on button press in pushbutton10.
function pushbutton10_Callback(hObject, eventdata, handles)  
% in dialog screen, atoms in initial crystal of handles.x will be drawed
% and the view point is along the at&el direction 
at = str2num(get(handles.edit91,'string'));
el = str2num(get(handles.edit92,'string'));
step1forshowcrystal(at, el, hObject, eventdata, handles);

function step1forshowcrystal(at, el, hObject, eventdata, handles) 
handles.x_new(:,2:4) = handles.x(:,2:4);       % read coordinates from initial crystal     
Drawsupercell(hObject, eventdata, handles)
axes(handles.axes1)
hold on;
view(at, el);  %view along the at&el direction


%the coordinates of atoms in projection will be edited, which will be
%sliced by multi-slice method
atompos=handles.x(:,2:4);  % read atoms' coordinators
% rotate atoms via RR matrix
at=-at;
el=90+el;

RR=[cosd(at) -sind(at) 0; cosd(el)*sind(at)  cosd(el)*cosd(at)  -sind(el);
sind(at)*sind(el) sind(el)*cosd(at) cosd(el)];
newatom=(RR*atompos.').';

%these coordinates will be sliced in the following codes 
atompos(:,3) = newatom(:,3)-min(newatom(:,3));  % the beginning height is 0 
atompos(:,1) = newatom(:,1)-min(newatom(:,1))+str2num(get(handles.edit100,'string'));
atompos(:,2) = newatom(:,2)-min(newatom(:,2))+str2num(get(handles.edit101,'string'));
handles.x_new(:,2:4)=atompos;
handles.x_mysavenewz = atompos(:,3);  %remember z to slicing crystal better ��z���ֵ��¼ס�������ڻ��ֲ�󣬿�����ʱ��ʾ��ȷ��ԭ�Ӹ߶ȣ�������ֲ㲻���⣬�´��ٷֲ�͵�������߶�

Drawsupercell_figure(handles.x_new, handles.dis_atomsize, handles.x_mysavenewz, 0, 0)
title('Rotate crystal and view from top(Red-X, Green-Y, Blue-Z)')
view([0 0 -1])  
guidata(hObject, handles);

% old codes about rotation of crystal are in v5 version.



% --- Executes on button press in pushbutton16.
function pushbutton16_Callback(hObject, eventdata, handles)
%according to the screen show to project crystal
[at,el] = view();
set(handles.edit91,'string',num2str(at));
set(handles.edit92,'string',num2str(el));
step1forshowcrystal(at, el, hObject, eventdata, handles);



function edit91_Callback(hObject, eventdata, handles)
function edit91_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



% --- Executes on selection change in popupmenu4.
function popupmenu4_Callback(hObject, eventdata, handles)
function popupmenu4_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function nx=mycen(x);
if rem(x,2)==0;
    nx=x/2+1;
else
    nx=(x+1)/2;
end



function edit94_Callback(hObject, eventdata, handles)
function edit94_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit1_Callback(hObject, eventdata, handles)
function edit1_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit99_Callback(hObject, eventdata, handles)
function edit99_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit98_Callback(hObject, eventdata, handles)
function edit98_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --------------------------------------------------------------------
function Untitled_23_Callback(hObject, eventdata, handles)

function Untitled_24_Callback(hObject, eventdata, handles)
prompt='Atom Size';
dlg_title='Atomsize( constant>0 ):';
num_lines = 1;%����Ի��������;
default_val={ num2str(handles.dis_atomsize) };%Ĭ�ϵ�ֵ;
answer=inputdlg(prompt,dlg_title,num_lines,default_val);
handles.dis_atomsize=str2num(answer{1});

guidata(hObject, handles);


% --------------------------------------------------------------------
function Untitled_25_Callback(hObject, eventdata, handles)
Drawsupercell(hObject, eventdata, handles)

function edit100_Callback(hObject, eventdata, handles)
function edit100_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit101_Callback(hObject, eventdata, handles)
function edit101_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --------------------------------------------------------------------
function Untitled_26_Callback(hObject, eventdata, handles)  %��ʾ�ֲ��ͼ��
cd(handles.tempcd)  %������һ�δ򿪵�Ŀ¼
[FileName,PathName]=uigetfile({'*.mat','Slicing Crystal File (*.mat)'});%���ļ�
if length(FileName)==1 & FileName(1)==0;  %���ѡ��cancalû���ļ��Ļ���������ֹ
    cd(handles.execd)
    return;
end
load(strcat(PathName, FileName));
cd(handles.execd)
%load allresul.mat
%������ά��slicing��ƽ�棻
Drawsupercell_figure(allresul.x_new, handles.dis_atomsize, allresul.x_new(:,4), 0, 0);  %
hold on;   

maxx = max(allresul.x_new(:,2));
maxy = max(allresul.x_new(:,3));
[yy,xx] = meshgrid(linspace(0, maxy, max(10, round( maxy/ (20* handles.dis_atomsize)))), linspace(0, maxx, max(10, round( maxx/ (20* handles.dis_atomsize)))));
for i=0:length( allresul.slicethick )-1
    mesh(xx,yy,i*allresul.eachthick.*ones(size(xx)));
end

%%�ٶ໭һ����棬��ɫ�ģ���Ϊ�����棻
surf(xx,yy,(i+1)*allresul.eachthick.*ones(size(xx)), zeros(size(xx)));

disp('Project on the top slice')
title('Sices of 3D show');


% --- Executes on button press in pushbutton13.
function pushbutton13_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_26 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% ypj
if exist('paraslicing.txt','file') == 0
    x=[10 0 0 0.2 1];
    save paraslicing.txt -ascii x
end
load paraslicing.txt
myvar.x_new = handles.x_new;
myvar.slicenum = paraslicing(1);  %�⼸����ť����Ҫ�ˡ�
myvar.top = paraslicing(2);
myvar.bottom = paraslicing(3);
myvar.x_mysavenewz = handles.x_mysavenewz;
disp(strcat('From top to bottom, distance is :', num2str(max(handles.x_mysavenewz)-min(handles.x_mysavenewz))))
%Kepeng
myvar.atomsize = paraslicing(4);
myvar.Projected = paraslicing(5);
% myvar.Pick = paraslicing(6);
%Kepeng
slicingdisplay(myvar);
guidata(hObject, handles);


% --- Executes on button press in radiobutton13.
function radiobutton13_Callback(hObject, eventdata, handles)
if get(handles.radiobutton13,'value')==1 
    set(handles.radiobutton10, 'value', 0);
    set(handles.radiobutton9, 'value', 0);
    set(handles.radiobutton11, 'value', 0);
end
if get(handles.radiobutton13,'value')==1
    set(handles.text115, 'visible', 'off');  %HRTEM��عر�
    
    set(handles.text21, 'visible', 'off');   % aper��mrad��1/nm��λ
    set(handles.text125, 'visible', 'on');
    set(handles.text24, 'visible', 'on');   % convergence & aperture��mrad��1/nm��λ
    set(handles.text134, 'visible', 'off');
    
    set(handles.text106, 'visible', 'on');  %STEM���
    set(handles.text108, 'visible','on');
    set(handles.text112, 'visible','on');
    set(handles.text110, 'visible','on');
    set(handles.edit77, 'visible','on');
    set(handles.edit79, 'visible','on');
    
    set(handles.text28, 'visible','on');
    set(handles.text103, 'visible','on');
    set(handles.edit20, 'visible','on');
    set(handles.edit19, 'visible','on');
        
else
    set(handles.text106, 'visible', 'off');
    set(handles.text108, 'visible','off');
    set(handles.text112, 'visible','off');
    set(handles.text110, 'visible','off');
    set(handles.edit77, 'visible','off');
    set(handles.edit79, 'visible','off');
    
    set(handles.text28, 'visible','off');
    set(handles.text103, 'visible','off');
    set(handles.edit20, 'visible','off');
    set(handles.edit19, 'visible','off');
    
    set(handles.text21, 'visible', 'on');   % aper��mrad��1/nm��λ
    set(handles.text125, 'visible', 'off');
    set(handles.text24, 'visible', 'off');   % convergence & aperture��mrad��1/nm��λ
    set(handles.text134, 'visible', 'on');
    
end
guidata(hObject, handles);



function edit105_Callback(hObject, eventdata, handles)
function edit105_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function radiobutton14_Callback(hObject, eventdata, handles)



function edit106_Callback(hObject, eventdata, handles)
disp(strcat('Results will be saved in the folder:', handles.saveresult));
disp(strcat('Name as:', get(handles.edit106,'string'),'_*'));
guidata(hObject, handles);

function edit106_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton18.
function pushbutton18_Callback(hObject, eventdata, handles)
selpath = uigetdir();  % select one dir to save results
handles.saveresult = selpath;
disp(strcat('Results will be saved in the folder:', handles.saveresult));
disp(strcat('Name as:', get(handles.edit106,'string'),'_*'));
guidata(hObject, handles);


% --------------------------------------------------------------------
function Untitled_27_Callback(hObject, eventdata, handles)

prompt = {'Source Size (Angstrom)'; 'Sampling (Integer number)'};
dlg_title = 'Converlution for STEM images';
num_lines = 1;
def = {num2str(handles.conv_source);num2str(handles.conv_sampling)};
answer = inputdlg(prompt,dlg_title,num_lines,def);
if isempty(answer) 
    return; 
end
handles.conv_source=str2num(answer{1}); 
handles.conv_sampling=str2num(answer{2}); %�Ѿ���ĵ�Դ�ߴ磬�Լ��Ŵ��ʳ���


cd(handles.saveresult)  %������һ�δ򿪵�Ŀ¼
[FileName,PathName]=uigetfile({'*.dat','DATA File (*.dat)';'*.*','All File (*.*)'});%���ļ�
if length(FileName)==1 & FileName(1)==0;  %���ѡ��cancalû���ļ��Ļ���������ֹ
    cd(handles.execd)
    return;
end

cd(handles.execd);  %�ص��ɵ��ļ���
sx = str2num(get(handles.edit19, 'string'));
sy = str2num(get(handles.edit20, 'string'));

fid = fopen(strcat(PathName, FileName), 'r');
tempimage = fread(fid, [sy, sx], 'float');
fclose(fid);
tempimage= tempimage.';

sampling = str2num(get(handles.edit75, 'string'));% image sampling rate
step = str2num(get(handles.edit77, 'string'));  %step

if handles.conv_sampling>0
    [oldx, oldy]=size(tempimage);
    newx = round(oldx*handles.conv_sampling);
    newy = round(oldy*handles.conv_sampling);
    newimage = zeros(newx, newy);
    tempfft = fftshift(fft2(tempimage))/oldx/oldy;
    newimage(mycen(newx)-mycen(oldx)+1:mycen(newx)-mycen(oldx)+oldx, mycen(newy)-mycen(oldy)+1:mycen(newy)-mycen(oldy)+oldy) = tempfft;

    if handles.conv_source>0  %the input image sampling rate is sampling*step; the zoomed image sampling rate is sampling*step/handles.conv_sampling
        [x,y] = meshgrid( ([1:newy]-mycen(newy))/(newy*sampling*step/handles.conv_sampling) ,([1:newx]-mycen(newx))/(newx*sampling*step/handles.conv_sampling));
        ss = handles.conv_source/(sqrt(-log(0.5))*2);
        newimage = newimage.*exp(-((pi*ss)^2*(x.^2+y.^2)));
    else
        disp('Input a value (source size) larger than 0')
        return;
    end

    newimageresul = real( ifft2(ifftshift(newimage))*newx*newy);
else
    disp('Input a value (converlution sampling) larger than 0')
    return;
end

figure;imshow(newimageresul, 'XData',...
           str2num(get(handles.edit17,'string')) + [1:newy-1].*sampling*step/handles.conv_sampling, ...
           'YData', str2num(get(handles.edit18,'string')) + [0:newx-1].*sampling*step/handles.conv_sampling,...
          'DisplayRange',[]);axis on;
      
tempname = FileName; tempname(strfind(tempname,'_')) = ' ';
title(strcat('File name:', tempname))

fid=fopen(strcat(PathName, FileName(1:end-4),'_source_',num2str( handles.conv_source), '_zoom',num2str(handles.conv_sampling),'.dat'), 'w');
fwrite(fid, newimageresul.', 'float');
fclose(fid);
figure;surf(newimageresul)
shading interp
view([0 -1 0])
guidata(hObject, handles);
% 
% NyO = round(handles.Oversampling*Ny);%Ny,NxΪͼ��ߴ�
% NxO = round(handles.Oversampling*Nx);
%     imgft = zeros(NyO,NxO);
%     NxOMid = floor(handles.Oversampling*Nx/2)+1;
%     NyOMid = floor(handles.Oversampling*Ny/2)+1;
%     NxMid = floor(Nx/2)+1;
%     NyMid = floor(Ny/2)+1;
%     imgft(NyOMid-NyMid+[1:Ny],NxOMid-NxMid+[1:Nx]) = fftshift(fft2(img));
%     % Apply effective source size for STEM images:
%     if ss > 0%ss��sorcesize�й�,ss=sourcesize     ss = ss/(sqrt(-log(0.5))*2);
%         [qx,qy] = meshgrid((-NxOMid+[1:NxO])/(Nx*dx),(-NyOMid+[1:NyO])/(Ny*dy));
%         img = ifft2(ifftshift(imgft.*exp(-((pi*ss)^2*(qx.^2+qy.^2)))));        
%     else
%         img = ifft2(ifftshift(imgft));
%     end
%     if realFlag
%        img = real(img); 
%     end
%     scale = 1/(Nx*Ny);
%     Nx = round(handles.Oversampling*Nx);
%     Ny = round(handles.Oversampling*Ny);
%     dx = dx/handles.Oversampling;
%     dy = dy/handles.Oversampling;
%     scale = (scale*Nx*Ny);
% else
%     if ss > 0
%         NxMid = floor(Nx/2)+1;
%         NyMid = floor(Ny/2)+1;
%         [qx,qy] = meshgrid((-NxMid+[1:Nx])/(Nx*dx),(-NyMid+[1:Ny])/(Ny*dy));
%         img = real(ifft2(fft2(img).*ifftshift(exp(-((pi*ss)^2*(qx.^2+qy.^2))))));        
%     end    
% end


% --- Executes on selection change in popupmenu5.
function popupmenu5_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu5 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu5
if get(handles.popupmenu5,'value')==3  %�������lobato��ϵ�����Ͳ��ֲܷ�
    set(handles.checkbox1, 'value', 0);
end
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function popupmenu5_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in CPURB.
function CPURB_Callback(hObject, eventdata, handles)
% hObject    handle to CPURB (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if get(handles.CPURB, 'value')
    set(handles.GPURB,'value',0);
else
    set(handles.GPURB,'value',1);
end
guidata(hObject, handles);

% Hint: get(hObject,'Value') returns toggle state of CPURB


% --- Executes on button press in GPURB.
function GPURB_Callback(hObject, eventdata, handles)
% hObject    handle to GPURB (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if get(handles.GPURB, 'value')
    set(handles.CPURB,'value',0);
else
    set(handles.CPURB,'value',1);
end
guidata(hObject, handles);
% Hint: get(hObject,'Value') returns toggle state of GPURB


% --- Executes during object creation, after setting all properties.
function pushbutton3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pushbutton3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --------------------------------------------------------------------
function Untitled_28_Callback(hObject, eventdata, handles)
cd(handles.tempcd)  %������һ�δ򿪵�Ŀ¼
[FileName,PathName]=uigetfile({'*.cif','CIF File (*.cif)'});%���ļ�
if length(FileName)==1 & FileName(1)==0;  %���ѡ��cancalû���ļ��Ļ���������ֹ
    cd(handles.execd)
    return;
end
set(handles.edit1,'string',FileName);

lenxyz = [0,  0,  0];
prompt={sprintf(strcat('Cube''s width heigh thickness',num2str(lenxyz)))};
dlg_title = 'Input';
num_lines = 1;
defaultans = {'0';'0';'0'};

lenxyz=inputdlg(prompt,dlg_title,num_lines,defaultans);
lenxyz = str2num( cell2mat(lenxyz) );       

        
    z=unique(handles.x(:,1));
    for i=1:length(z)
        prompt={sprintf('DW factor for the %.0f Atom ',z(i))};
        dlg_title = 'Input';
        num_lines = 1;
        defaultans = {'0'};
        DW=inputdlg(prompt,dlg_title,num_lines,defaultans);
        handles.x(find(handles.x(:,1)==z(i)),7)=str2num(cell2mat(DW));
    end
    temp = handles.x;
    xlswrite(strcat(PathName, FileName(1:end-4)),temp);
    disp(strcat('Parameters of atoms are saved in the file', strcat(PathName, FileName(1:end-4)),'.xlsx'));


cd(handles.execd);  %�ص��ɵ��ļ���
handles.tempcd=PathName;  %��¼��һ���򿪵��ļ���

%��ȡpdb����cif�ľ�������
if sum(FileName(end-2:end)=='xls')==3 || sum(FileName(end-3:end)=='xlsx')==4 %��ʾ��pdb��ʽ
    handles.x=xlsread(strcat(PathName, FileName));
else
    disp('Cannot read this file!');
    return;
end

atompos=handles.x(:,2:4);  % ��λ�䰣~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
atompos(:,3) = atompos(:,3)-min(atompos(:,3)); %+1.5*eachthick;  % ԭ�Ӹ߶ȣ���Сֵ��Ϊ0��������һ����ղ��ȣ�����ֲ�
                                                            % ���ұ�֤����λ�ڵ�0�㡣��Ϊ����������λ�ڵڶ��㣬������ȡ��floor��������
                                                            % ��������3,0-3-6-9������1.5*eachslice��ֵΪ4.5�����ڵ�һ�����м�

atompos(:,1) = atompos(:,1)-min(atompos(:,1)); %+handles.extendis;
atompos(:,2) = atompos(:,2)-min(atompos(:,2)); %+handles.extendis;

handles.atompos = atompos;  %�����޸�һ��ԭ�ӵ����꣬������Ϣ��Ҫ�鿴handles.x�����ȡ�Ľ����
                %ֻ��handles.x�е�������Ҫ�ı�һ�£�ֻ�������ò��ʱ��
handles.x(:,2:4)=atompos;  %���¸�ԭ�ӵ����꣬��֤���Ǵ���0�ģ����б߽�
%%%%%%%%%%%%%%%%%%%%%%%-----------------end 202012290930

% the atoms recoded in handles.x_new is just for draw. 
% handles.x_new will be reedit according to the view point and crystal
% rotation. 
%Also see the function of 'Rotation & view' and 'Projected along view'
%---begin 202012290935
handles.x_new=handles.x;  %����תǰ�Ľṹ��¼������handles.x_new��¼������ת��
Drawsupercell(hObject, eventdata, handles)
guidata(hObject, handles);


% --------------------------------------------------------------------
function Untitled_29_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_29 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function Untitled_30_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_30 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function Untitled_31_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_31 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function Untitled_32_Callback(hObject, eventdata, handles)
%���벨��������Ҫ�ߴ���Ϣ��
[FileName,handles.batchs.PathName]=uigetfile({'*.dat','Data Only Format (*.dat)'});%�򿪺�׺Ϊtxt���ļ�
if length(FileName)==1 & FileName(1)==0;  %���ѡ��cancalû���ļ��Ļ���������ֹ
    mycddir(handles.execd)
    return;
end
%ѡ�񲨺����󣬲�������ͼ��Ҫ��ʾ
    prompt={'Enter wave size: width (Cubic)';...
       'Input flag: Show wave(S) and Run simulation(R)'; 'Images'' number of ';'Process number';...
         'Batch size for one process';...
         'TOP-LEFT coordinate (x-y)';...  %��Ҫһ������������ֵ
         'Intercept size (width & hight)' };   %��Ҫһ������������ֵ
    name='Input Image Info.';
    numlines=1;
    defaultanswer={num2str(handles.batchs.wavesx),...
        handles.batchs.showorrun,...
        num2str(handles.batchs.totalnumber),...
        num2str(handles.batchs.processnum),...
        num2str(handles.batchs.GPUBatch),...
        num2str(handles.batchs.top_left),...
        num2str(handles.batchs.width_heigh)};
    answer=inputdlg(prompt,name,numlines,defaultanswer);
    options.Resize='on';
    options.WindowStyle='normal';

    handles.batchs.wavesx=str2num(answer{1});handles.batchs.wavesy=handles.batchs.wavesx;
    handles.batchs.showorrun = answer{2};
    handles.batchs.totalnumber = str2num(answer{3});
    handles.batchs.processnum = str2num(answer{4});
    handles.batchs.GPUBatch = str2num(answer{5});
    handles.batchs.top_left = str2num(answer{6});
    handles.batchs.width_heigh = str2num(answer{7});
    
   

fid=fopen(strcat(handles.batchs.PathName,FileName),'r');  %��ȡ����������ʯīϩ  %1204��sampling��Ϊ0.1������0.03������̫����ͼ���С��320
a=fread(fid,[handles.batchs.wavesx*2,handles.batchs.wavesy],'float');%
fclose(fid);
disimage=a(1:2:end,:)+sqrt(-1)*a(2:2:end,:);   %
mywave=disimage.';

if handles.batchs.showorrun == 'S'  %ֻ���ȿ�һ����ѡ�������Ƿ���ʡ�
    %�Ѳ�������ʾ���������һ�����ȡ������λ��
    figure;imshow(angle(disimage),[]); title('phase of wave');   %�Ѳ�������ʾ����������Ҫ��ȡ���������Ͻ��ú�Ȧ�����½��û��Ǳ�ʾ
    hold on; plot(handles.batchs.top_left(2),handles.batchs.top_left(1),'ro');   
    plot(handles.batchs.top_left(2)+handles.batchs.width_heigh(1),handles.batchs.top_left(1)+handles.batchs.width_heigh(2),'y*');
    %��ʾ�Խ���
    line([handles.batchs.top_left(2), handles.batchs.top_left(1)], [handles.batchs.top_left(2)+handles.batchs.width_heigh(1), handles.batchs.top_left(1)+handles.batchs.width_heigh(2)],'color','r')
    %��ʾ������
    axis on
    return;
end

%�洢ͼ����ļ�������
handles.batchs.savePathName=uigetdir();


%���ּ���Ͷ����ݿ�ʼ
sx=handles.batchs.wavesx;
sy=handles.batchs.wavesy;
sxy=handles.batchs.wavesx

%��ȡһЩ���������ѹ���������޹�
h_ba=(6.6255916e-34)/2/pi;
h=6.6255916e-34;
e=1.602102e-19;
me=9.1093807e-31;
c=2.9979251e+8;
vol=str2num(get(handles.edit_Vol,'string'));
handles.lambda=1.0e+9*h.*c./sqrt(e*vol*1000*(2*me*c*c+e*vol*1000));  %���㲨������λnm
handles.vol=vol*1000;
handles.sampling=str2num(get(handles.edit75,'string'));  %��/pixel  

%%��ɫ�ĳߴ磬��ɫ����ʸ��gx���ڼ���ԭ��λ�õģ�s2ʹ��peng�Ĳ�������ԭ���Ƴ��ֲ�
[gx_green,gy_green]=meshgrid((-sxy/2:(sxy/2-1))./(sxy*handles.sampling), ...
    (-sxy/2:(sxy/2-1))./(sxy*handles.sampling)); %REPRO��λ��1/A)^2,����(1/nm)^2
sx_green=gx_green/2;  % peng ��s����
sy_green=gy_green/2;
s2_green=sx_green.^2+sy_green.^2;

mywave=fftshift(fft2(mywave));   %���첨���������������ȡҲ����;��Ƶ�ʿռ�


%��ȡ��ɢ�Ĳ�����
[para_part1, U, V]=CommonPara_TEM(handles, sx, sy);  %��CTF��PhasePlate�ù̶��Ĵ�С����
U=U-para_part1.tiltx;  %add 20210119 the incident wave is tilted but not the CTF is tilted
V=V-para_part1.tilty;  %the tilted CTF is shown for CTF display. but in simulation, only a tilted incident beam is required.
para_part2 = readtccfromhandles_newTEM(hObject, handles, para_part1.lambda);   %ֻ��һ���������ѿ���е����ݶ���ȡ�ˣ�������֤�����ط���Ҫ��ȡ����ʱ����ͳһ�Ĵ���

myap=ones(size(U));
myap=myaperture(myap,1/(para_part1.sampling.*sx),1/(para_part1.sampling.*sy),0,para_part1.gmax,0,0,1);  
handles.gfsf=para_part1.gfsf;

%Ϊÿ����ֵ���÷�Χ
%��ȡA1-��������ֵ�ķ���%�Լ��Ƕȷ�Χ
[txtName,txtPathName]=uigetfile({'*.txt','TXT Only Format (*.txt)'});%�򿪺�׺Ϊtxt���ļ�
if length(txtName)==1 & txtName(1)==0;  %���ѡ��cancalû���ļ��Ļ���������ֹ
    mycddir(handles.execd)
    return;
else
    %��ȡ���в�������ֵ
%newpara_part2.phiA1=newpara_part2.phiA1+(rand-0.5)*2*myrandrange.phiA1;
%%�����ɣ����еĽǶȣ�Ҳ����һ����Χ�ڣ��Ӽ�һ����ֵ������Ӽ�180����������180��һȦ
allvalue = load(strcat(txtPathName,txtName));
myrandrange.focus=allvalue(1); %4.8;%7.5; %units: nm will be chanced with both positive and nagative directions
myrandrange.A1=allvalue(2);%0.96;%6; %units:nm
myrandrange.phiA1=allvalue(3);%0.96;%6; %units:nm
myrandrange.A2=allvalue(4);%10;%49.268;%30; %units:nm
myrandrange.phiA2=allvalue(5);%0.96;%6; %units:nm
myrandrange.B2=allvalue(6);%10;%15;%16.4;%15; %units:nm ****  �м�*�Ķ���0709�޸ĵķ�Χ
myrandrange.phiB2=allvalue(7);%0.96;%6; %units:nm
myrandrange.Cs=allvalue(8);%1000; %3000;%2220;%3000; %unit nm  ****
myrandrange.A3=allvalue(9);%700; %50;%2220;%50; %units:nm  ****
myrandrange.phiA3=allvalue(10);%0.96;%6; %units:nm
myrandrange.S3=allvalue(11);%600; %100;%558;%100; %units:nm ****
myrandrange.phiS3=allvalue(12);%allvalue(3);%0.96;%6; %units:nm
myrandrange.A4=allvalue(13);%50000; %5000;%94800;%5000; %units:nm ****
myrandrange.phiA4=allvalue(14);%0.96;%6; %units:nm
myrandrange.D4=allvalue(15);10000; %20000;%19000;%20000; %units:nm ****
myrandrange.phiD4=allvalue(16);%0.96;%6; %units:nm
myrandrange.B4=allvalue(17);%40000; %5000;%19000;%5000; %units:nm ****
myrandrange.phiB4=allvalue(18);%0.96;%6; %units:nm
myrandrange.A5=allvalue(19);%1000000; %5000000;%3800000;%5000000; %units:nm ****
myrandrange.phiA5=allvalue(20);%0.96;%6; %units:nm
myrandrange.C5=allvalue(21);%1000000; %2000000;%3800000;%2000000; %units:nm ****  %��ʵҲ�������ðٷֱȡ�

% myrandrange.focus=4;%4.8;%7.5; %units: nm will be chanced with both positive and nagative directions
% myrandrange.A1=2.5;%0.96;%6; %units:nm
% myrandrange.A2=10;%49.268;%30; %units:nm
% myrandrange.B2=10;%15;%16.4;%15; %units:nm ****  �м�*�Ķ���0709�޸ĵķ�Χ
% myrandrange.Cs=1000; %3000;%2220;%3000; %unit nm  ****
% myrandrange.A3=700; %50;%2220;%50; %units:nm  ****
% myrandrange.S3=600; %100;%558;%100; %units:nm ****
% myrandrange.A4=50000; %5000;%94800;%5000; %units:nm ****
% myrandrange.D4=10000; %20000;%19000;%20000; %units:nm ****
% myrandrange.B4=40000; %5000;%19000;%5000; %units:nm ****
% myrandrange.A5=1000000; %5000000;%3800000;%5000000; %units:nm ****
% myrandrange.C5=1000000; %2000000;%3800000;%2000000; %units:nm ****  %��ʵҲ�������ðٷֱȡ�

% for k=1:91
%      A1_phi(k)=-184+4*k;  
% end
% A2_phi=A1_phi;
% B2_phi=A1_phi;
% 
% for k=1:37
%      A3_phi(k)=-190+10*k;
% end
% S3_phi=A3_phi;
% A4_phi=A3_phi;
% B4_phi=A3_phi;
% D4_phi=A3_phi;
% A5_phi=A3_phi;

end
%rand_value=0*rand_value;

yita=str2num(get(handles.edit_Spr,'String')); 


%������������ı仯
beginnumold=0;
tic

    
myap=ones(size(gx_green));
myap=myaperture(myap,1/(handles.sampling.*sxy),1/(handles.sampling.*sxy),0,0.1*str2num(get(handles.edit_Ape,'string')),0,0,1);
mywave=mywave.*myap;

commpara=[para_part1.sampling, vol, ... %�����ʣ�kv��λ�ĵ�ѹ
        para_part1.gmax, yita, para_part1.mtotal, para_part1.alafa*1000, ...  %���Ƶ�ʣ�focus spread������ĸ�˹�����ĸ����� beam convergence
        para_part1.tilt*1000, -para_part1.phitilt];
    
rand('seed',sum(100*clock));

othervar = [handles.batchs.processnum, handles.batchs.GPUBatch,...
    handles.batchs.totalnumber, handles.batchs.top_left, handles.batchs.width_heigh];
savePathName = handles.batchs.savePathName;
save('pydata.mat', 'para_part1', 'para_part2', 'myrandrange', ...
    'U', 'V', 'mywave', 'commpara', 'othervar', 'savePathName');
dos("python cuda_stem64.py");

% neixunhuan=6;
% for j=1:100
%     beginnum=beginnumold+(j-1)*neixunhuan;
%     %allimage = zeros(28*44,neixunhuan);
%     %allpara = zeros(33, neixunhuan);
%     
%   %  rand_value=(rand(1000,21)-0.5)*2;  %���Ʋ����ĸ���������2��һ���
%     for i=1:neixunhuan; %length(rand_value(:,1))
%         newpara_part2=para_part2;
%         %newpara_part2.focus=round(newpara_part2.focus+(rand-0.5)*2*myrandrange.focus,1);  %��������������ұ���С�������λ
%         newpara_part2.focus=newpara_part2.focus+(rand-0.5)*2*myrandrange.focus; %add at 20200820��ʵ����
%         %newpara_part2.A1=round(newpara_part2.A1+(rand-0.5)*2*myrandrange.A1,1);  %������
%         newpara_part2.A1=newpara_part2.A1+(rand-0.5)*2*myrandrange.A1; %add at 20200820��ʵ����
%         newpara_part2.phiA1=A1_phi(round(rand*(90))+1); %ʵ����
%         %newpara_part2.phiA1=round(rand*180); %add at 20200820��������
%         %newpara_part2.A2=round(newpara_part2.A2+(rand-0.5)*2*myrandrange.A2,0);  %������
%         newpara_part2.A2=newpara_part2.A2+(rand-0.5)*2*myrandrange.A2; %add at 20200820��ʵ����
%         newpara_part2.phiA2=A2_phi(round(rand*(90))+1); %ʵ����
%         %newpara_part2.phiA2=round(rand*180); %add at 20200820
%         %newpara_part2.B2=round(newpara_part2.B2+(rand-0.5)*2*myrandrange.B2,0);  %��������������ұ���С�����1λ
%         newpara_part2.B2=newpara_part2.B2+(rand-0.5)*2*myrandrange.B2; %add at 20200820,ʵ����
%         newpara_part2.phiB2=B2_phi(round(rand*(90))+1); %ʵ����
%         %newpara_part2.phiB2=round(rand*180); %add at 20200820��������
%         %newpara_part2.Cs=round(newpara_part2.Cs+(rand-0.5)*2*myrandrange.Cs,-2);  %��������������ұ���С�����1λ
%         newpara_part2.Cs=newpara_part2.Cs+(rand-0.5)*2*myrandrange.Cs; %add at 20200820��ʵ����
%         %newpara_part2.A3=round(newpara_part2.A3+(rand-0.5)*2*myrandrange.A3,-2);  %��������������ұ���С�����1λ
%         newpara_part2.A3=newpara_part2.A3+(rand-0.5)*2*myrandrange.A3; %add at 20200820��ʵ����
%         newpara_part2.phiA3=A3_phi(round(rand*(36))+1); %ʵ����
%         %newpara_part2.phiA3=round(rand*180); %add at 20200820��������
%         %newpara_part2.S3=round(newpara_part2.S3+(rand-0.5)*2*myrandrange.S3,-2);  %��������������ұ���С�����1λ
%         newpara_part2.S3=newpara_part2.S3+(rand-0.5)*2*myrandrange.S3; %add at 20200820��ʵ����
%         newpara_part2.phiS3=S3_phi(round(rand*(36))+1); %ʵ����
%         %newpara_part2.phiS3=round(rand*180); %add at 20200820��������
% %����Ĳ������Ͳ������        %newpara_part2.A4=round(newpara_part2.A4+(rand-0.5)*2*myrandrange.A4,-3);  %��������������ұ���С�����1λ
%         %newpara_part2.phiA4=A4_phi(round(rand*(36))+1);
%         %newpara_part2.D4=round(newpara_part2.D4+(rand-0.5)*2*myrandrange.D4,-3);  %��������������ұ���С�����1λ
%         %newpara_part2.phiD4=D4_phi(round(rand*(36))+1);
%         %newpara_part2.B4=round(newpara_part2.B4+(rand-0.5)*2*myrandrange.B4,-3);  %��������������ұ���С����ǰ4λ
%         newpara_part2.B4=newpara_part2.B4+(rand-0.5)*2*myrandrange.B4; %add at 20200820��ʵ����
%         newpara_part2.phiB4=B4_phi(round(rand*(36))+1); %ʵ����
%         %newpara_part2.phiB4=round(rand*180); %������
%        %newpara_part2.A5=round(newpara_part2.A5+(rand-0.5)*2*myrandrange.A5,-5);  %��������������ұ���С����ǰ4λ
%       %newpara_part2.phiA5=A5_phi(round(rand*(36))+1);
%         %newpara_part2.C5=round(newpara_part2.C5+(rand-0.5)*2*myrandrange.C5,-5);  %��������������ұ���С����ǰ4λ
% 
%     
%       newpara_part2.A1x=newpara_part2.A1*cos(newpara_part2.phiA1/180*pi);   newpara_part2.A1y=newpara_part2.A1*sin(newpara_part2.phiA1/180*pi);
%     %C1=focus;  focus���������defocs������ƫ���
%       newpara_part2.A2x=newpara_part2.A2*cos(newpara_part2.phiA2/180*pi);   newpara_part2.A2y=newpara_part2.A2*sin(newpara_part2.phiA2/180*pi);
%       newpara_part2.B2x=newpara_part2.B2*cos(newpara_part2.phiB2/180*pi);   newpara_part2.B2y=newpara_part2.B2*sin(newpara_part2.phiB2/180*pi);
%         newpara_part2.A3x=newpara_part2.A3*cos(newpara_part2.phiA3/180*pi);   newpara_part2.A3y=newpara_part2.A3*sin(newpara_part2.phiA3/180*pi);
%        newpara_part2.S3x=newpara_part2.S3*cos(newpara_part2.phiS3/180*pi);   newpara_part2.S3y=newpara_part2.S3*sin(newpara_part2.phiS3/180*pi);
%        newpara_part2.C3=newpara_part2.Cs;
%       newpara_part2.A4x=newpara_part2.A4*cos(newpara_part2.phiA4/180*pi);   newpara_part2.A4y=newpara_part2.A4*sin(newpara_part2.phiA4/180*pi);
%       newpara_part2.D4x=newpara_part2.D4*cos(newpara_part2.phiD4/180*pi);   newpara_part2.D4y=newpara_part2.D4*sin(newpara_part2.phiD4/180*pi);
%       newpara_part2.B4x=newpara_part2.B4*cos(newpara_part2.phiB4/180*pi);   newpara_part2.B4y=newpara_part2.B4*sin(newpara_part2.phiB4/180*pi);
%       newpara_part2.A5x=newpara_part2.A5*cos(newpara_part2.phiA5/180*pi);   newpara_part2.A5y=newpara_part2.A5*sin(newpara_part2.phiA5/180*pi);
%       newpara_part2.S5x=newpara_part2.S5*cos(newpara_part2.phiS5/180*pi);   newpara_part2.S5y=newpara_part2.S5*sin(newpara_part2.phiS5/180*pi);
%        newpara_part2.C5=newpara_part2.C5;
%        newpara_part2.D5x=newpara_part2.D5*cos(newpara_part2.phiD5/180*pi);   newpara_part2.D5y=newpara_part2.D5*sin(newpara_part2.phiD5/180*pi);
%    
%        % para_part2 = readtccfromhandles_newTEM(hObject, handles, para_part1.lambda);%����������ϵ����ݣ�������һ��
%        % newpara_part2=para_part2; %����������ϵ����ݣ�������һ��
%        
%       %  residual_aberr=WholeTCC_residual_newTEM(newpara_part2, U,V);  %����phase plate visualing the residual aberrations
%         mytcc=WholeTCC2D_newTEM_forpar(para_part1, newpara_part2, U,V);
%     
% 
%        %��ĳ�ʼ��
%        myinten=zeros(sxy,sxy);
%        %����������
% %        for nn=1:para_part1.mtotal;
% %           myinten=myinten+handles.gfsf(nn).*abs(ifft2(ifftshift(mywave.*mytcc(:,:,nn).*myap))).^2;
% %        end
%        for nn=1:para_part1.mtotal;
%           myinten=myinten+handles.gfsf(nn).*abs(ifft2(ifftshift(mywave.*mytcc(:,:,nn)))).^2;  %��myap����ֱ�Ӱ�myap����ȥ��
%        end
%        
%        savimage=zeros(28,44,3);
%        myresul=myinten(124:151,145:188); 
%        myresul=(myresul-min(myresul(:)))./(max(myresul(:))-min(myresul(:)));
%        savimage(:,:,1)=myresul;
%        savimage(:,:,2)=myresul;
%        savimage(:,:,3)=myresul;
%        %20201208����
%        imwrite(savimage, strcat(num2str(beginnum+i),'.jpg'));
%        
%        %allimage(:,i)= myresul(:);
%    
%  
% 
%     
% %     %�������sin
% %     %myresul = sin(myresul*255./10);  %%%%%%%%%%%%%%%%0619
% %     myresul = sin(myresul*25);  %%%%%%%%%%%%%%%%0619
% %     myresul = (myresul-min(myresul(:)))./(max(myresul(:))-min(myresul(:))); %%%%%%%%%%%%%%%%0619
% %    
% %     savimage=repmat(myresul,1,1,3);
% %     %savimage(:,:,1)=myresul;savimage(:,:,2)=myresul;savimage(:,:,3)=myresul;
% %     imwrite(savimage,strcat('C:\Users\EM_Lab\Desktop\cos_0820_2i\',num2str(beginnum+i),'.jpg')); 
%     
%     myallpara=[commpara, ...
%         newpara_part2.focus, ...
%         newpara_part2.A1, -newpara_part2.phiA1, ...
%         newpara_part2.A2, -newpara_part2.phiA2, ...
%         newpara_part2.B2, -newpara_part2.phiB2, ...
%         newpara_part2.Cs/1000, ...
%         newpara_part2.A3/1000, -newpara_part2.phiA3, ...
%         newpara_part2.S3/1000, -newpara_part2.phiS3, ...
%         newpara_part2.A4/1000, -newpara_part2.phiA4, ...
%         newpara_part2.D4/1000, -newpara_part2.phiD4, ...        
%         newpara_part2.B4/1000, -newpara_part2.phiB4, ...
%         newpara_part2.C5/1000000, ...
%         newpara_part2.A5/1000000, -newpara_part2.phiA5, ...
%         newpara_part2.D5/1000000, -newpara_part2.phiD5, ...        
%         newpara_part2.S5/1000000, -newpara_part2.phiS5];  %ͨ��tccread����Ĳ������Ƕ���Ҫȡ��
%         
%         %allpara(:,i) = myallpara;
%        
%         fid=fopen(strcat(num2str(beginnum+i),'.para'),'w');
%         fwrite(fid, myallpara, 'float');
%         fclose(fid);
%     end
%         %fid=fopen(strcat('H:\goku\1204\batch3\',num2str(j),'img.dat'),'w');
%         %fwrite(fid, allimage, 'float');
%         %fclose(fid);    
%         %fid=fopen(strcat('H:\goku\1204\batch3\',num2str(j),'.para'),'w');
%         %fwrite(fid, allpara, 'float');
%         %fclose(fid);
% end


%figure;imshow(myinten,[]);  %��ʾͼ��
toc

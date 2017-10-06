function varargout = mbc_wmm(varargin)
% MBC_WMM MATLAB code for display_vmfmm_samples.fig
%      MBC_WMM, by itself, creates a new MBC_WMM or raises the existing
%      singleton*.
%
%      H = MBC_WMM returns the handle to a new MBC_WMM or the handle to
%      the existing singleton*.
%
%      MBC_WMM('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in MBC_WMM.M with the given input arguments.
%
%      MBC_WMM('Property','Value',...) creates a new MBC_WMM or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before mbc_wmm_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to mbc_wmm_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help mbc_wmm

% Last Modified by GUIDE v2.5 15-Sep-2014 14:41:09

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @mbc_wmm_OpeningFcn, ...
                   'gui_OutputFcn',  @mbc_wmm_OutputFcn, ...
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


% --- Executes just before mbc_wmm is made visible.
function mbc_wmm_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to mbc_wmm (see VARARGIN)

% Choose default command line output for mbc_wmm
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes mbc_wmm wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = mbc_wmm_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in load_data_GT.
function load_data_GT_Callback(hObject, eventdata, handles)
% hObject    handle to load_data_GT (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Select the appropriate file
[FileName,PathName,FilterIndex] = uigetfile;

% Load samples and labels
load(strcat(PathName, FileName));

% global data 
handles.wmmSample = wmmSample;
handles.labels = labels;
guidata(hObject,handles);

% Display data in sphere with label
showWatsonSamples_gui(handles, wmmSample, labels, length(unique(labels))/2);

% Enable clustering button
set(handles.cluster_data,'Enable','on');

% --- Executes on key press with focus on load_data_GT and none of its controls.
function load_data_GT_KeyPressFcn(hObject, eventdata, handles)
% hObject    handle to load_data_GT (see GCBO)
% eventdata  structure with the following fields (see UICONTROL)
%	Key: name of the key that was pressed, in lower case
%	Character: character interpretation of the key(s) that was pressed
%	Modifier: name(s) of the modifier key(s) (i.e., control, shift) pressed
% handles    structure with handles and user data (see GUIDATA)


% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over load_data_GT.
function load_data_GT_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to load_data_GT (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in load_data_only.
function load_data_only_Callback(hObject, eventdata, handles)
% hObject    handle to load_data_only (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Select the appropriate file
[FileName,PathName,FilterIndex] = uigetfile;

% Load samples
load(strcat(PathName, FileName));

% global data 
handles.wmmSample = wmmSample;
guidata(hObject,handles);

% Display data in sphere without label
showWatsonSamples_gui(handles, wmmSample, ones(1,size(wmmSample,1)), 1);

% Enable clustering button
set(handles.cluster_data,'Enable','on');

% --- Executes on button press in load_data_GT.
function load_data_GT_Callback_Callback(hObject, eventdata, handles)
% hObject    handle to load_data_GT (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in load_data_only.
function load_data_only_Callback_Callback(hObject, eventdata, handles)
% hObject    handle to load_data_only (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in save_plot.
function save_plot_Callback(hObject, eventdata, handles)
% hObject    handle to save_plot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Select file name
FileName = uiputfile('*.*');
export_fig(handles.axes1, strcat(FileName, '.png'));


% --- Executes during object creation, after setting all properties.
function figure1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
addpath('exportfig');
addpath('vmfmatlab');


% --- Executes on button press in gen_vmfmm.
function gen_vmfmm_Callback(hObject, eventdata, handles)
% hObject    handle to gen_vmfmm (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

cla(handles.axes1,'reset');

% Generate samples from the mixture model
numSamp = 10000;
num_comp = str2num(get(handles.tb_numclass, 'String'));
[wmmSample, labels, ~] = Samples_Watson_mixture_model_gui(num_comp, numSamp);

% global data 
handles.wmmSample = wmmSample;
handles.labels = labels;
guidata(hObject,handles);

% Display data in sphere without label
showWatsonSamples_gui(handles, wmmSample, labels, num_comp);

% Enable clustering button
set(handles.cluster_data,'Enable','on');

% --- Executes on button press in cluster_data.
function cluster_data_Callback(hObject, eventdata, handles)
% hObject    handle to cluster_data (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Load necessary data from global data 
handles = guidata(hObject);
wmmSample = handles.wmmSample;

% Perform clustering data
kMax = 15;

% get model parameters with k_max
[~, params, ~, ~] = getWMMParams_bd(wmmSample, kMax, 'dm');

% remove invalid indices
params = removeInvalidIndices(params);

% get model parameters with k_max
[allIC, allClust, ~] = wmm_HC_BD_IC(wmmSample, params);

% Model selection
y = allIC.BIC;
x = 1: length(y);

% [~, numComp] = min(y); % BIC
% [~, numComp] = min(clInfo.icInfo.Beta_min); % Beta min
% [~, numComp] = min(clInfo.icInfo.ICL); % ICL
[numComp,~, ~] = l_method(x,y); % L-method
% [numComp,~, ~] = reg_pcws2(x,y,[1 30],[1 1]); % WLR-30

% Get the final clustering results
finalClust = allClust(:, numComp);

% Display data in sphere with label
spread_gui(handles.axes1, wmmSample', finalClust);

% Disable clustering button
set(handles.cluster_data,'Enable','off');



function tb_numclass_Callback(hObject, eventdata, handles)
% hObject    handle to tb_numclass (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of tb_numclass as text
%        str2double(get(hObject,'String')) returns contents of tb_numclass as a double


% --- Executes during object creation, after setting all properties.
function tb_numclass_CreateFcn(hObject, eventdata, handles)
% hObject    handle to tb_numclass (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

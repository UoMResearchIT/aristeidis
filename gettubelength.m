function out = gettubelength
    f = figure('Units','Normalized',...
        'Position',[.4 .4 .3 .3],...
        'NumberTitle','off',...
        'Name','Info');
    set(gcf,'Color',[0.6,0.3,0.8])
    e = uicontrol('Style','Edit',...
        'Units','Normalized',...
        'Position',[.1 .4 .3 .1],...
        'Tag','myedit');
    p1= uicontrol('Style','PushButton',...
        'Units','Normalized',...
        'Position',[.6 .4 .3 .1],...
        'String','Done',...
        'CallBack','uiresume(gcbf)');
    p2= uicontrol('Style','text',...
        'Units','Normalized',...
        'Position',[.3 .2 .4 .05],...
        'String','Give length of heating tube (in m) ',...
        'CallBack','uiresume(gcbf)');
    uiwait(f)
    out = str2num(get(e,'String'));
    close(f)
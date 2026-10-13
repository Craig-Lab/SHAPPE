function IQR = plotIQR(data,figHandle,color)

    figure(figHandle)
    IQR = quantile(data(~isnan(data)), [0.25 0.75]);
    hStat=fill([0,length(data)+1,length(data)+1,0],[IQR(1),IQR(1),IQR(2),IQR(2)],color,'EdgeColor','none');
    hStat.FaceAlpha=0.75;

end
function  mydraw(Pop,prombeg,promend)%画个体图，总体图
%Pop的第一条即为sequence
    f=[];
    hold off;
    delete('G:/iGEM/IOfile/Rfile/*')
    path='G:/iGEM/IOfile/proms.xlsx';
    delete(path);
    fastaWrite(Pop,path);
    system('"G:\R\R\R-4.1.0\bin\R.exe" --no-save --slave < G:\iGEM\R\fitness2.R');
    affinities=load('G:/iGEM/IOfile/Rfile/output.txt');
    learningcurve=affinities(1:25,1);
    for i=1:size(affinities,2) % 重置评分
        [f(i),affinities(:,i)]=seqarea_0(learningcurve,affinities(:,i),prombeg-73,promend-73);
    end
    %% 个体图
    standardx=[0 1113];
    srandardy=[0 0];
    y=affinities;
    x=[1:length(y)];
    hold on
    for i=1:size(y,2)
        plot(x,y(:,i),standardx,srandardy,'-r','LineWidth',1);
    end
    axis([0 1113 -14 4]);
    xlabel('bpnum');
    ylabel('affinity');
    title('Traditional-parameter=nan');
   
    
    
%     for i=1:size(affinities,2) % 自动画图并保存
%         y=affinities(:,i);
%         x=[1:length(y)];
%         plot(x,y,'-b',standardx,srandardy,'-r','LineWidth',4);
%         xlabel('bpnum');
%         ylabel('affinity');
%         title(num2str(i));
%         grid on;
%         set(gca,'GridLineStyle','--','GridColor','k','GridAlpha',0.5,'LineWidth',2);
%         axis([0 1113 -20 4]);
%         saveas(gcf,['G:/iGEM/IOfile/Rfile/',num2str(i),'.jpg']); 
%     end
    %% 总体图
%     y=f;
%     x=[1:length(y)];
%     plot(x,y,'-b','LineWidth',4);
%     grid on;
%     set(gca,'GridLineStyle','--','GridColor','k','GridAlpha',0.5,'LineWidth',2);
%     axis([0 length(f) -10000 10]);
%     title('Overall afinity');
%     xlabel('generation')
%     ylabel('affinity')
%     saveas(gcf,['G:/iGEM/IOfile/Rfile/','0verall','.jpg']);
end


function [resiinfo, series_resi] = ...
             CalResidualInfo(all_n, all_n_i, eachthick, slicethick, multiceng);  
         
         series_n=[]; series_n_i=[]; %ÿ������˶��ٵ�ԭ��
ele_n=[]; ele_n_i=[]; 
newabsorp_n=[]; newabsorp_n_i=[];


if ~isempty(all_n)
    for i=1:length(slicethick)
        jj = find(all_n(:,4) >= slicethick(i)-multiceng*eachthick & all_n(:,4) < slicethick(i)+multiceng*eachthick+eachthick);  %�ҵ��������Ҫ���ԭ�����
        series_n(i) = length(jj);  %��¼�ж��ٸ�ԭ��λ����һ����
        ele_n(end+1:end+length(jj),1:16) = all_n(jj,:);
        
         %��������10.bʽ��
        z1(1:length(jj),1) = ele_n(end-series_n(i)+1:end, 4) -(slicethick(i)+eachthick);
        z2(1:length(jj),1) = ele_n(end-series_n(i)+1:end, 4) - slicethick(i);
        if i==1
            z2(:)=inf;
        end
        if i==length(slicethick)
            z1(:)=-inf;  
        end
        z2(find(abs(z2)>=multiceng*eachthick)) = inf;
        z1(find(abs(z1)>=multiceng*eachthick)) = -inf;
        ele_n(end-series_n(i)+1:end, 17) = project_poten_contribute(ele_n(end-series_n(i)+1:end,8),  z1, z2);
        ele_n(end-series_n(i)+1:end, 18) = project_poten_contribute(ele_n(end-series_n(i)+1:end,10),  z1, z2);
        ele_n(end-series_n(i)+1:end, 19) = project_poten_contribute(ele_n(end-series_n(i)+1:end,12),  z1, z2);
        ele_n(end-series_n(i)+1:end, 20) = project_poten_contribute(ele_n(end-series_n(i)+1:end,14),  z1, z2);
        ele_n(end-series_n(i)+1:end, 21) = project_poten_contribute(ele_n(end-series_n(i)+1:end,16),  z1, z2);
        
        
        if ~isempty(absorp_n)  %����Ƿ�Ҫ��������
            newabsorp_n(end+1:end+length(jj),1:10) = absorp_n(jj,:);
            newabsorp_n(end-series_n(i)+1:end, 11) = project_poten_contribute(newabsorp_n(end-series_n(i)+1:end,2),  z1, z2);
            newabsorp_n(end-series_n(i)+1:end, 12) = project_poten_contribute(newabsorp_n(end-series_n(i)+1:end,4),  z1, z2);
            newabsorp_n(end-series_n(i)+1:end, 13) = project_poten_contribute(newabsorp_n(end-series_n(i)+1:end,6),  z1, z2);
            newabsorp_n(end-series_n(i)+1:end, 14) = project_poten_contribute(newabsorp_n(end-series_n(i)+1:end,8),  z1, z2);
            newabsorp_n(end-series_n(i)+1:end, 15) = project_poten_contribute(newabsorp_n(end-series_n(i)+1:end,10),  z1, z2);
        end
        z1=[];  %��Ҫ���㣬�������ݿ��ܻᴮ
        z2=[];
    end
    %��֮�󣬰���Ҫ�����
    ele_n(:,7)=ele_n(:,7).*ele_n(:,17);
    ele_n(:,9)=ele_n(:,9).*ele_n(:,18);
    ele_n(:,11)=ele_n(:,11).*ele_n(:,19);
    ele_n(:,13)=ele_n(:,13).*ele_n(:,20);
    ele_n(:,15)=ele_n(:,15).*ele_n(:,21);
    ele_n(:,17:21) = [];

    ele_n(:,5)=[]; %ɾ��DB
    ele_n(:,4)=[]; %ɾ��Z
    ele_n(:,1)=[]; %ɾ��ԭ������
    
    %����֮�ϣ�������Ϊ��
    % ele_n�У��ֱ��ǣ�
    % ԭ������������1-3��DW���ӣ�ԭ����*ռ���ʣ�ab��������ţ�ÿ����˹�ı���
    %֮�����Ϊ������������
    % ele_n������13��
    %����1-2��ռ����*����3��abϵ��   
 
    if ~isempty(absorp_n)  %����Ƿ�Ҫ��������
        %��֮�󣬰���Ҫ�����
       newabsorp_n(:,1)=newabsorp_n(:,1).*newabsorp_n(:,11);
       newabsorp_n(:,3)=newabsorp_n(:,3).*newabsorp_n(:,12);
       newabsorp_n(:,5)=newabsorp_n(:,5).*newabsorp_n(:,13);
       newabsorp_n(:,7)=newabsorp_n(:,7).*newabsorp_n(:,14);
       newabsorp_n(:,9)=newabsorp_n(:,9).*newabsorp_n(:,15);

       newabsorp_n(:,11:15)=[];
    end
end
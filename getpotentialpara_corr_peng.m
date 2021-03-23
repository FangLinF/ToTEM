function [corr_nuclear,corr_nuc_ion]  = getpotentialpara_corr_peng(x);  %如果需要导入修正势场的参数
load ree_l12p5.mat
re_neu=ree_l12p5;   %把质子性的写进去

%把系数导入到程序里，第一列是原子序号，2-3是a系数，4-5是b系数，newss在这个程序没有，是拟合的起点；
%newstart和newlast是Q函数的起点和终点s分量
% re=[xu, aaar, bbbr, newss,newstart,newlast,a1_newf,a1,a1_newf_peng,a2_newf,a2,a2_newf_peng,...
%     a3_newf,a3,a3_newf_peng,a4_newf,a4,a4_newf_peng,a5_newf,a5,a5_newf_peng,a6_newf,a6,a6_newf_peng];
re_ion=[];  %现在没有ion的离子性  建议第一列和第二列分别是原子序号和价态，后面类推一个位置

corr_nuclear=[]; corr_nuc_ion=[];  %结果先赋空
corr_ion=[];

%case 1，只有原子性或只有离子性
b=find(x(:,6)==0);  %找出只有原子性的原子
if ~isempty(b)
    corr_nuclear(:,1) = x(b,1);  %获得质子数
    corr_nuclear(:,2:4) = x(b,2:4); %获得坐标
    corr_nuclear(:,5) = x(b,7); %DB
    corr_nuclear(:,6) = x(b,8).*(1-x(b,6)); %占据率*离子性
    corr_nuclear(:,7:12) =  re_neu(corr_nuclear(:,1), [2:5,7:8]); %每个原子的修正系数导入
end

b=find(x(:,6)==1);  %找出只有离子性的原子
len=length(b);  %个数有几个离子性的  %the ion parameters comes from peng's
if ~isempty(b)
    corr_ion(:,1) = x(b,1);  %获得质子数
    corr_ion(:,2) = x(b,5);  %获得价态
    corr_ion(:,3:5) = x(b,2:4); %获得坐标
    corr_ion(:,6) = x(b,7); %DB
    corr_ion(:,7) = x(b,8).*x(b,6); %占据率*离子性
    
    [tempjiatai, aa, xuhao] =unique(re_ion(:,1:2),'rows','stable');
    clear aa
    %把质子数和价态对应的参数，都赋值到新的参数里面
    for i=1:length(tempjiatai(:,1))
        k = find(tempjiatai(i,1)==re_ion(:,1) & tempjiatai(i,2)==re_ion(:,2));
        tempionpara(i, 8:13)  = re_ion(k,[3:6,8:9]);  %
    end
    corr_ion(:,8:17) = tempionpara(xuhao, :);   %把所有原子的离子性都赋值。从tempionpara中提取
end

%这里if的判断，得到了单纯离子性或原子性的规定矩阵。详见 程序说明.doc文件
if ~isempty(corr_ion)   
    corr_nuclear(end+1:end+len,:)=corr_ion(:,[1,3:end]);
    clear corr_ion
end

% case 2
%找出既有原子性又有离子性的原子。---------------------------------
b=find(x(:,6)>0 & x(:,6)<1);  %找出只有离子性的原子
if ~isempty(b)
    corr_nuc_ion(:,1) = x(b,1);  %获得质子数
    corr_nuc_ion(:,2) = x(b,5);  %获得价态
    corr_nuc_ion(:,3:5) = x(b,2:4); %获得坐标
    corr_nuc_ion(:,6) = x(b,7); %DB
    corr_nuc_ion(:,7) = x(b,8); %占据率
    corr_nuc_ion(:,8) = x(b,6); %离子性
    
    [tempjiatai, aa, xuhao] =unique(corr_nuc_ion(:,1:2),'rows','stable');
    clear aa
    %把质子数和价态对应的参数，都赋值到新的参数里面
    for i=1:length(tempjiatai(:,1))
        k = find(tempjiatai(i,1)==re_ion(:,1) & tempjiatai(i,2)==re_ion(:,2));
        tempionpara(i, 1:6)  = re_ion(k,[3:6,8:9]);
    end
    corr_nuc_ion(:,9:14) = tempionpara(xuhao, :);   %把所有原子的离子性都赋值。从tempionpara中提取
    
    corr_nuc_ion(:,15) = 1- x(b,6); %原子性
    corr_nuc_ion(:,16:21) = re_neu(corr_nuc_ion(:,1), [2:5,7:8]);  %每个原子的散射参数
    corr_nuc_ion(:,2)=[];  %总共只剩下20列
end
return;
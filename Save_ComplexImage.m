function Save_ComplexImage(myI,filename)
[Ncol,Nrow]=size(myI);
%���渴��ͼ��ķ�����
temp=myI.';  %�洢��������ķ���  ��Ҫͼ����ת��
tempall(1:2:2*Nrow-1,:)=real(temp);
tempall(2:2:2*Nrow,:)=imag(temp);
fid=fopen(filename,'w');
fwrite(fid,tempall,'float');
fclose(fid);
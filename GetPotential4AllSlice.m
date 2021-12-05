   function potential=GetPotential4AllSlice(green_Ncol, green_Nrow,... 
    ele_n, absorp_n, ....�������ӻ���ԭ�ӵĵ��Ժ�����
    ele_n_i, absorp_n_i, ... %ԭ�ӻ�ԭ��+���ӣ����Ի�������
    series_n, series_n_i, ...  %ԭ�����д���
    s2, gx_green, gy_green, ...
    sigma, PARAMETER, APERTURE);   %Ϊ��STEM���㲻����������ֻ���뵽HRTEM��CBED��
potential=zeros(green_Ncol, green_Nrow, max(length(series_n), length(series_n_i)));  %ÿ���Ƴ���������
k1=0; %��¼���ĸ�ԭ�ӿ�ʼ�����Ƴ���
k2=0;
   %����ele_n�ľ���
if ~isempty(series_n)
    for i=1:length(series_n);  %����ÿһ����Ƴ�    
        for atom_i=k1+1:sum(series_n(1:i));
            k1=sum(series_n(1:i));  %Ϊ��һ�㣬���¼����ԭ�ӵ���ʼ����
            if ~isempty(ele_n)&isempty(absorp_n)
               %��������ԭ�ӻ����ӵĵ����Ƴ�����û�п��Ƿֲ����;
                potential(:,:,i)=potential(:,:,i)+PARAMETER*ele_n(atom_i,3).*...
                   ( ele_n(atom_i,4).*exp(-s2.*ele_n(atom_i,5)) ...
                   +ele_n(atom_i,6).*exp(-s2.*ele_n(atom_i,7)) ...
                   +ele_n(atom_i,8).*exp(-s2.*ele_n(atom_i,9)) ...
                   +ele_n(atom_i,10).*exp(-s2.*ele_n(atom_i,11)) ...
                   +ele_n(atom_i,12).*exp(-s2.*ele_n(atom_i,13)) ) ...
                   .*exp(-2*pi*sqrt(-1)*(gx_green*ele_n(atom_i,1) + gy_green*ele_n(atom_i,2)));
            elseif ~isempty(ele_n)&~isempty(absorp_n)
                 potential(:,:,i)=potential(:,:,i)+PARAMETER*ele_n(atom_i,3)*...
                  ( ele_n(atom_i,4).*exp(-s2.*ele_n(atom_i,5)) ...
                   +ele_n(atom_i,6).*exp(-s2.*ele_n(atom_i,7)) ...
                   +ele_n(atom_i,8).*exp(-s2.*ele_n(atom_i,9)) ...
                   +ele_n(atom_i,10).*exp(-s2.*ele_n(atom_i,11)) ...
                   +ele_n(atom_i,12).*exp(-s2.*ele_n(atom_i,13))  ...
                   +sqrt(-1)*( absorp_n(atom_i,1).*exp(-s2.*absorp_n(atom_i,2)) ...
                   +absorp_n(atom_i,3).*exp(-s2.*absorp_n(atom_i,4)) ...
                   +absorp_n(atom_i,5).*exp(-s2.*absorp_n(atom_i,6)) ...
                   +absorp_n(atom_i,7).*exp(-s2.*absorp_n(atom_i,8)) ...
                   +absorp_n(atom_i,9).*exp(-s2.*absorp_n(atom_i,10)) )) ...
                   .*exp(-2*pi*sqrt(-1)*(gx_green*ele_n(atom_i,1) + gy_green*ele_n(atom_i,2)));
            end
        end
    end
end
if ~isempty(series_n_i)
    for i=1:length(series_n_i);  %����ÿһ����Ƴ�    %֮ǰ©�������ѭ��
      for atom_i=k2+1:sum(series_n_i(1:i));
          k2=sum(series_n_i(1:i));  %Ϊ��һ�㣬���¼����ԭ�ӵ���ʼ����
          if ~isempty(ele_n_i) & isempty(absorp_n_i)    %����ԭ�Ӻ������ԵĹ��ױ���
              potential(:,:,i)=potential(:,:,i)+PARAMETER* ele_n_i(atom_i,3)*(ele_n_i(atom_i,4).*...
                   ( ele_n_i(atom_i,5).*exp(-s2.*ele_n_i(atom_i,6)) ...
                   +ele_n_i(atom_i,7).*exp(-s2.*ele_n_i(atom_i,8)) ...
                   +ele_n_i(atom_i,9).*exp(-s2.*ele_n_i(atom_i,10)) ...
                   +ele_n_i(atom_i,11).*exp(-s2.*ele_n_i(atom_i,12)) ...
                   +ele_n_i(atom_i,13).*exp(-s2.*ele_n_i(atom_i,14)) ) ...
                   + ele_n_i(atom_i,15).*...
                   ( ele_n_i(atom_i,16).*exp(-s2.*ele_n_i(atom_i,17)) ...
                   +ele_n_i(atom_i,18).*exp(-s2.*ele_n_i(atom_i,19)) ...
                   +ele_n_i(atom_i,20).*exp(-s2.*ele_n_i(atom_i,21)) ...
                   +ele_n_i(atom_i,22).*exp(-s2.*ele_n_i(atom_i,23)) ...
                   +ele_n_i(atom_i,24).*exp(-s2.*ele_n_i(atom_i,25)) )) ...
                   .*exp(-2*pi*sqrt(-1)*(gx_green*ele_n_i(atom_i,1) + gy_green*ele_n_i(atom_i,2)));
          elseif ~isempty(ele_n_i) & ~isempty(absorp_n_i)
                   potential(:,:,i)=potential(:,:,i)+PARAMETER* (ele_n_i(atom_i,3)*ele_n_i(atom_i,4).*...
                   ( ele_n_i(atom_i,5).*exp(-s2.*ele_n_i(atom_i,6)) ...
                   +ele_n_i(atom_i,7).*exp(-s2.*ele_n_i(atom_i,8)) ...
                   +ele_n_i(atom_i,9).*exp(-s2.*ele_n_i(atom_i,10)) ...
                   +ele_n_i(atom_i,11).*exp(-s2.*ele_n_i(atom_i,12)) ...
                   +ele_n_i(atom_i,13).*exp(-s2.*ele_n_i(atom_i,14)) ) ...
                   + ele_n_i(atom_i,3)*ele_n_i(atom_i,15).*...
                   ( ele_n_i(atom_i,16).*exp(-s2.*ele_n_i(atom_i,17)) ...
                   +ele_n_i(atom_i,18).*exp(-s2.*ele_n_i(atom_i,19)) ...
                   +ele_n_i(atom_i,20).*exp(-s2.*ele_n_i(atom_i,21)) ...
                   +ele_n_i(atom_i,22).*exp(-s2.*ele_n_i(atom_i,23)) ...
                   +ele_n_i(atom_i,24).*exp(-s2.*ele_n_i(atom_i,25)) ) ...
                   + sqrt(-1)*( absorp_n_i(atom_i,1).*exp(-s2.*absorp_n_i(atom_i,2)) ...
                   +absorp_n_i(atom_i,3).*exp(-s2.*absorp_n_i(atom_i,4)) ...
                   +absorp_n_i(atom_i,5).*exp(-s2.*absorp_n_i(atom_i,6)) ...
                   +absorp_n_i(atom_i,7).*exp(-s2.*absorp_n_i(atom_i,8)) ...
                   +absorp_n_i(atom_i,9).*exp(-s2.*absorp_n_i(atom_i,10))))...
                   .*exp(-2*pi*sqrt(-1)*(gx_green*ele_n_i(atom_i,1) + gy_green*ele_n_i(atom_i,2)));
        end
      end
    end
end

for i=1:max(length(series_n), length(series_n_i))   %��Ҫ���໥���ó���
    potentialsave(:,:,i)=ifft2(ifftshift(potential(:,:,i).*APERTURE))*green_Ncol*green_Nrow;
    potential(:,:,i)=exp(sqrt(-1)*potentialsave(:,:,i)*sigma./1000);  %����kirkland��ֵ�Ĺ�ϵ��sigma�ĵ�λ��/ǧ��*A
end
return;

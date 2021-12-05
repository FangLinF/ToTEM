%cuda的核心计算
function myresul=cuda_STEM_core_complex(PARAMETER, sigma, ...%包含了哪些原子，每层原子多少个；参数parameter算potential的;相互作用系数sigama，算透过函数需要的
    ele_n, absorp_n, ....仅有离子或者原子的弹性和吸收
    ele_n_i, absorp_n_i, ... %原子或原子+离子，弹性或者吸收
    series_n, series_n_i, ...  %原子排列次序
    gx_green, gy_green, s2, APERTURE,...   %绿色区域计算potential时候，gx和gy是用来算原子位置相关的倒格矢，s2是算peng potential的,APERTURE是防止2/3光阑外发生wrap，详见笔记
    probe, gfsf, psf_fft, ...%probe的形式，还有几个矩阵的比例分配，传播矩阵
    green_Ncol, green_Nrow, ...%绿色的尺寸，绿色倒空间的sx和sy，
    probe_ingreenNrow, probe_ingreenNcol, ...
    width_red, hight_red, probestep, ...  %选取的相对位置的左上角x和y坐标，以及右下角x和y坐标
    APER);  %最后增加的光阑函数
potential=zeros(green_Ncol, green_Nrow, max(length(series_n), length(series_n_i)));  %每层势场都存下来
k1=0; %记录从哪个原子开始计算势场。
k2=0;
   %先算ele_n的矩阵
if ~isempty(series_n)
    for i=1:length(series_n);  %计算每一层的势场    
        for atom_i=k1+1:sum(series_n(1:i));
            k1=sum(series_n(1:i));  %为下一层，更新计算的原子的起始序数
            if ~isempty(ele_n)&isempty(absorp_n)
               %仅仅计算原子或离子的弹性势场，并没有考虑分层比例;
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
      for atom_i=k2+1:sum(series_n_i(1:i));
          k2=sum(series_n_i(1:i));  %为下一层，更新计算的原子的起始序数
          if ~isempty(ele_n_i) & isempty(absorp_n_i)    %考虑原子和离子性的贡献比例
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

for i=1:max(length(series_n), length(series_n_i))   %需要乘相互作用常数
    potential(:,:,i)=ifft2(ifftshift(potential(:,:,i).*APERTURE))*green_Ncol*green_Nrow;
    potential(:,:,i)=exp(sqrt(-1)*potential(:,:,i)*sigma./1000);  %根据kirkland数值的关系，sigma的单位是/千伏*A
end

%传播
[probesx, probesy]=size(probe);

myresul=zeros(hight_red, width_red, length(APER(1,1,:)));
for i=0:width_red-1;
    i
    for j=0:hight_red-1;
        for nn=1:length(gfsf);
            mywave=probe(:,:,nn);  %读入波函数
            for k=1:length(potential(1,1,:));
                mywave =  mywave.*potential(j*probestep+probe_ingreenNrow:j*probestep+probe_ingreenNrow+probesx-1, ...
                   i*probestep+probe_ingreenNcol:i*probestep+probe_ingreenNcol+probesx-1,k);
                   mywave=fft2(mywave);   %注意，节省掉一个fftshift
                   mywave=mywave.*psf_fft;
                   mywave=ifft2(mywave);
            end
            mywave=fftshift(fft2(mywave));
            for kk=1:length(APER(1,1,:))              
                myresul(j+1,i+1,kk)=myresul(j+1,i+1,kk)+gfsf(nn).*sum(sum(abs(mywave.*APER(:,:,kk)).^2));  %低通、带通，或高通函数
            end
        end
    end
end
% for kk=1:length(APER(1,1,:)) 
%     figure;imshow(myresul(:,:,kk)/handles.width_red/handles.height_red,[]);colorbar
% end
return
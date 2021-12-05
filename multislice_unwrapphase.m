function multislice_unwrapphase(mywave, potential, psf_fft, sampling, lambda, eachthick)
[e,f]=find(fftshift(abs(ifft2(psf_fft))) == max(max(fftshift(abs(ifft2(psf_fft))))))
%�ҵ�psf_fft������λ�á�
psf_realspace = fftshift(ifft2(psf_fft));
sizexy = 10;  %�涨һ��psf�ĳߴ�
smallpsf = psf_realspace(e-sizexy:e+sizexy, f-sizexy:f+sizexy);  %ȡһС���ֽ��о�� 
%������smallpsfֻ����Ϊ�ο�����Ϊwrap̫����ˣ����ѽ��
[sx,sy]=size(psf_fft); 
psf=zeros(sx, sy);  
%�����psf�����԰ɡ�������һ��������ʱ��psf��ͬ�ߴ�ľ�������ʱ��ʽ������ʱ��ĸ������󣬿����ܲ��ܶԵ��ϡ�
[xx,yy]=meshgrid([1:sx]-e, [1:sy]-f);
xx=xx*sampling;
yy=yy*sampling;
%����x��y�ľ���,R=(x,y),���¸���6.90ʽ��д��psf��ʱ��ı�ʾ
psf=(1/sqrt(-1)/lambda/eachthick)*exp(sqrt(-1)*pi/lambda/eachthick*(xx.^2+yy.^2));
%��ô�������psf������psf_realspace��һ������Ϊ�����ɳ������������psf_fft(find(p2>1/(16*handles.sampling*handles.sampling)))=0;
%���ԱȽϵ�ʱ�򣬿��ǰ������Ǿ䲻Ҫ�ˣ�Ȼ����бȽ�
%��ʵ��ͼ����ֱ��дpsf������Ҫ�����˲�����Ϊ�ڸ���Ҷ�������棬��������������ƵĻ�����ߺ��ұߣ��ϱߺ��±ߵģ�����ֻ��Ч�������������ͼ����ֱ��д����Ļ��������������ĸ��ţ���Χ��0�˵ġ�
%psf��Ӧ����������λ����ʾ��д��
%--------------------�õ�psfA��psfP
%����ʡ����һ�䣬��Ҫ�Ժ����д�İɡ�Ҫ��֤��psfsmall�ĸ���ֵһ��������λ���ڷ�ת


mywaveA = abs(mywave);
mywaveP = angle(mywave);
mywave = mywaveA.*exp(sqrt(-1)*mywaveP);  %�ڼ����Ƴ���ʱ��exp(sqrt(-1)*V)��ɢ���ƣ�
%����ԭ�����ģ������ƴ���λ����������ﶨ�壬�������е���λ���ֵ�������ֵ��������λ���
for k=1:length(potential(1,1,:));
        mywave =  mywave.*potential(:,:,k);
        
        potentialA = abs(potential(:,:,k));
        potentialP = unwrap( angle(potential(:,:,k)) );  
        %�Ѳ���������λ��ȷ�Ľ⿪����ʵ��Ӧ��GetPotential4AllSlice_multicore_lotabo_peng_corr
        %246�е�V���֡�ע�⣬Vֻ������ֵ�ġ�
        mywaveA = getmultiple_AP(mywaveA, mywaveP, potentialA, potentialP);   %�����������д
        %��mywaveA��mywaveP����Ϊ�˷���������Ϣ����ʾһ���������������������λ�����������Ƴ����
        
        [mywaveA, mywaveP] = getconv_AP(mywaveA, mywaveP, psfA, psfP);  %������˸����ĳ˷������и����ļӷ�����Ҫ���
%��ʽ����ͼ�������洦��ģ����Բ����ڸ���Ҷ�任�ˡ�
%         mywave=fft2(mywave);   %ע�⣬��ʡ��һ��fftshift
%         mywave=mywave.*psf_fft;
%         mywave=ifft2(mywave);
end
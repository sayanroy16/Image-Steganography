
clc;
clear all;
close all;

%code snippet for Character to Bit pattern conversion


message = 'hiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiii';
N_msg=length(message); %Character Length for Secret Message
N_int=8; %Size of the unsigned integer representation, let it be 8-bit
m1=dec2bin(message,N_int); % Converting the character to Unsigned integer form; Matrix of size N_msg x N_int
bn=zeros(1,N_msg*N_int); %Variable to store binary bit stream initialized with all bits 0.
m2=reshape(m1',1,N_msg*N_int); %Converting the matrix to a vector
pos=find((m2=='1')); %Logic to get index positions where bit value is 1 in message
bn(pos)=1;%Final bit stream to be processed further in next steps
 
 
 
%Extracting the R-plane from the image

c=imread('football.bmp'); %Reading RGB image
figure
imshow(c); %Displaying the RGB image
title('RGB Image');
c_r=c(:,:,1); %Here, c is an RGB image, means there are three planes R, G and B. To extract the planes, we use the index number : 1 indicates R, 2 indicates G and 3 indicates B. We extract 
%Red plane (R). Why R and why not G and B is to be checked??? Is there any
%change in we take G or B. Check it...
figure
imshow(c_r); %Displaying the R plane. Observe that it is a Gray Scale Image. Intensity values are 0,1,2,...255.
title('R Plane Image (Cover Image)');
[count,gl]=imhist(c_r); %Here we determine histogram for the gray scale image : gl gives the intensity values i.e.
% 0,1,2....255 and count gives the number of pixels in the image occupying
% a particular gray level. Ex. For gray value 0, the count of pixels is
% 24446 (Check it by typing in command window)

figure
imhist(c_r); %Plotting the histogram i.e. count vs gl (No. of pixels vs Gray Levels)
title('Histogram Plot of R-plane gray image');
xlabel('Gray Levels');
ylabel('No. of pixels');
 
%%%Next step proceed to extract Binary Image from R-plane gray intensity
%%%image. The input to LSB algorithm will be a. Binary Image b. Binary bit
%%%stream generated
 
%Gray R-plane image converted to binary image
 
%%My way (with thresholding)
 
th_ost=graythresh(c_r); %Determining the threshold for binary conversion of the image. We are using here OSTU's Global Threshold Technique for grayscale images.
%th_ost value will be between 0 and 1 always.
c_bin=imbinarize(c_r,th_ost); %Converting the grayscale image R-plane to binary image using OSTU thresholding.
 
 
figure
imshow(c_bin); %Observe the figure that the characteristics of the image are PRESERVED during binarization
title('Binary Image (My way)');
 
%Your way (No thresholding) 
 
% cd=double(c_r); % storing image information in cd 
% c1=mod(cd,2); % extracting single bit from the R plane
% 
% figure
% imshow(c1);%Observe the figure that the characteristics of the image are NOT PRESERVED during binarization
% title('Binary Image (Your way)');
 
%%For the next step when you go for LSB, you can check with both c_bin and
%%c1 images and see the results for each. I do not discard the way you have
%%implemented but feel that c_bin will perform better as there is
%%thresholding applied.
 
 
 
 
% %new part
% %encoding using LSB
% len=N_msg*N_int;%Each character will take 8 bits so total number of bits required to hide a single character is 8.  
% %binary_all='';  %binary_all will have the entire sequence of bits of the message
% [row,column]=size(c_bin); %returns the size of matrix c_r in separate variables column  and row.
% %We do this so that we can start two loops that runs through each of it?s elements and perform LSB.
% %stego=c_r;% we store c_r in a new variable stego, which will be the output image 
% stego=c_bin;
% % for i=1:length(message)% we apply a for loop from position 1 to the length of the message
% %     binary_all=[binary_all,m1(i,:)];% here we concatenate the entire sequence of message bits and m1,
% %     %which is a Matrix of size N_msg x N_int containing the binary form of message.
% % end
%  
% count=1;    %initializing count with 1
%  
% for i=1:row %we implement a for loop that runs from position 1 to end of the row
%     for j=1:column %we implement another for loop so that the program also runs through each element of all the rows
%         
%         %for every character in the message 
%         if count<=len 
%             %Obtain the LSB of the grey level of the pixel so as to compare the LSB with the message bit.
%             %When we use the modulo function what we do is get the remainder as the output. 
%             %So, when we divide 0 by 2 we get remainder as 0 and when we divide 1 by 2 we get the remainder as 1.
%             %Thus, we obtain the LSB which we will use to encode the message bits.
%             LSB=mod(c_bin(i,j),2);
%             
%             %Convert the bit from the message to numeric form so that we can perform bitwise XOR operation
% %             a=str2double(binary_all(count));
%  
%             %If the message bit and the LSB of the pixel are same,
%             %we set temp = 0 If the message bit and the LSB of the pixel are different, 
%             %we set temp = 1 This setting of temp can be done by taking XOR of message bit and the LSB of the pixel
%  
%             %temp=double(xor(LSB,a));
%             temp=xor(LSB,bn(count));
%             
%             %Update the pixel of output image to input image pixel value + temp and using the for loop we keep 
%             %updating the output image till all the bits in the message are embedded
%             stego(i,j)=c_bin(i,j)+temp;  
%             count=count+1; % we increase the count by 1
%         end
%     end
% end
%  
% 
%  
% %subplot(1,2,2);
% figure
% imshow(stego);% displaying the stego image (Binary)
% title('Stego Image (Binary Form)');
% 
% k=double(stego);
% 
% figure
% imshow(k);% displaying the stego image (Gray) 
% title('Stego Image (Gray Form)');

%LSB Encoding (follow this part)
len=N_msg*N_int; %Length of bit pattern to encode
stego=c_r; %Initialization : output image = image image
[r,c]=size(stego); %Rows and Columns of the image
count=1; % To count how many bits are encoded
for i=1:r
    for j=1:c
        if count<=len
            LSB=mod(stego(i,j),2);
            temp=xor(LSB,bn(count));
            stego(i,j)=stego(i,j)+uint8(temp);
            count=count+1; 
        end
    end
end

figure
imshow(stego);
title('Stego Image');

%Decoding the message
count=1;
message_in_bits='';
for i=1:r
    for j=1:c
        %For all the characters in the message
        if count<=len
            
            %Retrieve the LSB of the intensity level of the pixel
            LSB=mod(stego(i,j),2);
            
            %Append into message_in_bits to get bit sequence of message
            message_in_bits=[message_in_bits,num2str(LSB)];
        
            count=count+1;
        end
    end
end

%Converting the bit sequence into the original message
i=1;
original_message='';
while i<=len
    %Take a set of 8 bits at a time
    %Convert the set of bits to a decimal number
    %Convert the decimal number which is the ascii value to its corresponding character
    %Append the obtained character into the resultant string
    original_message=[original_message,char(bin2dec(message_in_bits(1,i:i+7)))];
    i=i+8;
end

disp(['The original message is: ',original_message]);


im1 = double(c_r); %initialization, gray image = im1
im2 = double(stego);%initialization, stego image = im2
%Find the mean squared error
mse = sum((im1(:)-im2(:)).^2) / prod(size(im1)); %the formula for calculating  mean squared error
% now find the psnr, as peak=255
psnr = 10*log10(255*255/mse); %the formula for calculating peak signal to noise ratio
 
val=ssim(stego,c_r);

figure
subplot(2,1,1)
stem(c_r,'r')
xlabel('Range of Pixels');
ylabel('No of Pixels');
title('Histogram of Original Image');
grid on

subplot(2,1,2)
stem(stego,'r')
xlabel('Range of Pixels');
ylabel('No of Pixels');
title('Histogram of stego Image');
grid on



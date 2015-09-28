% %region test
clear all;close all;
img = Get_file();
BW = im2bw(img);
[B,L,N,A] = bwboundaries(BW);
for k=1:length(B),
       if(~sum(A(k,:)))
         boundary = B{k};
         plot(boundary(:,2), boundary(:,1), 'r','LineWidth',2);
         for l=find(A(:,k))'
           boundary = B{l};
           plot(boundary(:,2), boundary(:,1), 'g','LineWidth',2);
         end
       end
end
figure;
imshow(BW); hold on;
     for k=1:length(B),
       boundary = B{k};
       if(k > N)
         plot(boundary(:,2), boundary(:,1), 'g','LineWidth',2);
       else
         plot(boundary(:,2), boundary(:,1), 'r','LineWidth',2);
       end
     end


%%
% seg = chenvese(img,'whole',800,0.2,'multiphase'); 
% 
% P = imread('anti-mass.jpg');
% % Imnoise the original input
% I = P;
% I(:,:,1) = imnoise(I(:,:,1),'speckle');
% I(:,:,2) = imnoise(I(:,:,2),'salt & pepper',0.8);
% figure(),subplot(1,2,1),imshow(P),title('original image');
% subplot(1,2,2),imshow(I),title('original image with two components adding noise')
% 
% % Normal Chan & Vese cannot work
% seg = chenvese(I,'large',300,0.02,'chan'); 
% 
% % Chan & Vese for vector image works here
% seg = chenvese(I,'large',300,0.02,'vector');
% % Using built-in mask = 'whole' leads faster and better segmentation
% seg = chenvese(I,'whole',800,0.02,'vector');
% % 
% % Customerlized Mask
% m = zeros(size(I,1),size(I,2));
% m(20:120,20:120) = 1;
% seg = chenvese(I,m,500,0.1,'chan'); % ability on gray image
% % Built-in Mask
% seg = chenvese(I,'medium',400,0.02,'chan'); % ability on gray image
% %-- End 
% 



%%
% [B,L,N,A] = bwboundaries(BW,'noholes');   
% figure,imshow(BW);hold on;
% 
% for k = 1:N
%     % Boundary k is the parent of a hole if the k-th column
%     % of the adjacency matrix A contains a non-zero element
%     if (nnz(A(:,k)) > 0)
%         boundary = B{k};
%         plot(boundary(:,2),...
%             boundary(:,1),'r','LineWidth',2);
%         % Loop through the children of boundary k
%         for l = find(A(:,k))'
%             boundary = B{l};
%             plot(boundary(:,2),boundary(:,1),'g','LineWidth',2);
%         end
%     end
% end

% 
% % 
% %  [B,L,N,A] = bwboundaries(BW);
% %      imshow(BW); hold on;
% %      colors=['b' 'g' 'r' 'c' 'm' 'y'];
% %      for k=1:length(B),
% %        boundary = B{k};
% %        cidx = mod(k,length(colors))+1;
% %        plot(boundary(:,2), boundary(:,1), colors(cidx),'LineWidth',2);
% %        %randomize text position for better visibility
% %        rndRow = ceil(length(boundary)/(mod(rand*k,7)+1));
% %        col = boundary(rndRow,2); row = boundary(rndRow,1);
% %        h = text(col+1, row-1, num2str(L(row,col)));
% %        set(h,'Color',colors(cidx),'FontSize',14,'FontWeight','bold');
% %      end
% %      figure; spy(A);
% try
% I = rgb2gray(img);    
% catch ex
%     I = im2bw(img);
% end
% %% get X Y
% ff = edge(I,'canny');
% [y, x] = find(ff == 1);
% edgesnap = [x(:),y(:)];
% set(0,'RecursionLimit',3000);

% a = sortedPoint(edgesnap);
% %[a,~] = dpsimplify(a,1);
% edgesnap = a;

% stat = regionprops(I,'BoundingBox');
% % plot(stat(:,1), stat(:,2));
% 
% 
% imcontour

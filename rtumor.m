function deltar = rtumor(coords,centerM)
 
  deltar=0;
  
  xpos(:,1)=coords(1,:);
  ypos(:,1)=coords(2,:);
  zpos(:,1)=coords(3,:);

 DT=delaunayTriangulation(xpos,ypos,zpos);

allF = [ DT(:,1) DT(:,2) DT(:,3); ...
    DT(:,1) DT(:,3) DT(:,4); ...
    DT(:,1) DT(:,4) DT(:,2); ...
    DT(:,2) DT(:,4) DT(:,3)];
 
 sortedF = sort(allF,2);

[u,m,n] = unique(sortedF,'rows');

counts = accumarray(n(:), 1);

 sorted_exteriorF = u(counts == 1,:);
 
 F = allF(ismember(sortedF,sorted_exteriorF,'rows'),:);
 
b = unique(F(:));

 for i=1:size(b,1)
   k=b(i,1);
    boundary(i,1)=xpos(k,1);
      boundary(i,2)=ypos(k,1);
        boundary(i,3)=zpos(k,1);
 end
 
 nBound=size(boundary,1);
 
 for k=1:nBound

       deltar= deltar + sqrt((boundary(k,1)-centerM(1,1))^2 + ...
       (boundary(k,2)-centerM(2,1))^2 + (boundary(k,3)-centerM(3,1))^2)/nBound;
   
end
end

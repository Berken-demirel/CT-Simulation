function y = read_txt_file()
[file,~,path] = uiputfile('*.txt','Workspace File');
fileID = fopen(file,'r');
A = fscanf(fileID,'%f');
n = size(A);
n = n(1,1);

a = A(1,1);
b = A(2,1);

projection_matrix = zeros(a,b);
for k = 0:(a-1)
    projection_matrix((k+1),:) = A((4+k * (b+1)) : (2+(k+1) * (b+1)));
end
 y = projection_matrix;
end
function newData = newData(X, N)

[sz1, sz2] = size(X);

classesSz = 24;

table = zeros(classesSz, sz2);
classes = zeros(classesSz+1, sz2);


for i=1:sz2
    t = hist(X(:,i),classesSz)/sz1;
    sz = (max(X(:,i)) - min(X(:,i)))/classesSz;
    classes(:,i) = min(X(:,i)):sz:max(X(:,i));
    table(:,i) = cumsum(t);
end

newData = zeros(N, sz2);

for i=1:N

num = rand(1);

for j=1:sz2
    value = find(table(:,j)>=num,1,"first");
    value = value(1);
    inter = [classes(value,j)
            classes(value+1,j)];
   newData(i,j) = unifrnd(inter(1), inter(2));
end
end
table
classes
end


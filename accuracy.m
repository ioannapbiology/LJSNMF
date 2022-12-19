function a = accuracy(x,y)
%compute the clustering evaluation metric 'Accuracy' of two vector
%clusterings
%   (sum of "hits" between the ground truth labels x and the mapped result
%   labels y)/# of data points. 
%   The best mapping can be found with the hungarian method


mapped_labels = bestMap(x, y);
%error_cnt = sum(x ~= partition_labels);

%get the Accuracy 'a'
a = length(find(x == mapped_labels))/length(x);
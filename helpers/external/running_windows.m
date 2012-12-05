function   [windowed_data]=running_windows(data,window_size, step_size)

%This function gets a vector of data and returns a matrix of
%N X (window_size). data is a 1 X N vector. The step size between the
%windows is step_size
%[windowed_data]=running_windows(data,window_size, step_size)

%        Written by Sigal Saar August 08 2005

size_data=size(data,2);

start_window=[1:floor(step_size):(size_data-window_size+1)]';
end_window=[window_size:floor(step_size):size_data]';

%adding zeros or sample point to make compensation for uneven number of lines watch out from doing a transpose to the matrix!!!
if (step_size-floor(step_size))~=0
    %step_size=floor(step_size);
    
    add_to_step=round(1/(step_size-floor(step_size)));
    
    n_data=zeros(1,length(start_window));
    
    
    add_to_window=[0:((size_data)/add_to_step)]'*ones(1,add_to_step+1);
    
    add_to_window=add_to_window';
    add_to_window=add_to_window(1:end)';
    %add_to_window=[0 ; add_to_window];
    
    start_window=start_window+add_to_window(1:length(start_window));
    end_window=end_window+add_to_window(1:length(start_window));
    
    
    good_index=find(end_window<size_data);
    
    start_window=start_window(good_index);
    end_window=end_window(good_index);
   

    n_data=n_data(1:length(end_window));
   n_data( add_to_step+1:add_to_step+1:end)=data(end_window(add_to_step+1:add_to_step+1:end)+1);
   if n_data(end)>size_data
       n_data(end)=0;
   end
end

char_string_dot=char(':'*ones(length(end_window),1));
char_string_coma=char(';'*ones(length(end_window),1));
starting_string=['[' ; char(' '*ones(length(end_window)-1,1))];
ending_string=[char(' '*ones(length(end_window)-1,1)) ; ']' ];

eval_matrix=[starting_string num2str(start_window) char_string_dot num2str(end_window) char_string_coma ending_string]';
if 0
 %Due to memory problems I had to split the calculations 
           matrix_is_1=eval([ eval_matrix(1:size(eval_matrix,1)*floor(size(eval_matrix,2)/8)) ']' ]);
           matrix_is_2=eval([ '[' eval_matrix(size(eval_matrix,1)*(floor(size(eval_matrix,2)/8)):size(eval_matrix,1)*(floor(size(eval_matrix,2)/8))*2) ']' ]);
           matrix_is_3=eval([ '[' eval_matrix(size(eval_matrix,1)*(floor(size(eval_matrix,2)/8))*2:size(eval_matrix,1)*(floor(size(eval_matrix,2)/8))*3) ']'  ]);
           matrix_is_4=eval([ '[' eval_matrix(size(eval_matrix,1)*(floor(size(eval_matrix,2)/8))*3:size(eval_matrix,1)*(floor(size(eval_matrix,2)/8))*4) ']' ]);
           matrix_is_5=eval([ '[' eval_matrix(size(eval_matrix,1)*(floor(size(eval_matrix,2)/8))*4:size(eval_matrix,1)*(floor(size(eval_matrix,2)/8))*5) ']' ]);
           matrix_is_6=eval([ '[' eval_matrix(size(eval_matrix,1)*floor(size(eval_matrix,2)/8)*5:size(eval_matrix,1)*(floor(size(eval_matrix,2)/8))*6) ']' ]);
           matrix_is_7=eval([ '[' eval_matrix(size(eval_matrix,1)*floor(size(eval_matrix,2)/8)*6:size(eval_matrix,1)*(floor(size(eval_matrix,2)/8))*7) ']' ]);
           matrix_is_8=eval([ '[' eval_matrix(size(eval_matrix,1)*(floor(size(eval_matrix,2)/8))*7:end)  ]);
   %        windowed_matrix=eval(tmp_char(1:(size(tmp_char,1)*size(tmp_char,2))));

windowed_data=data([ matrix_is_1 ; matrix_is_2 ;matrix_is_3 ;matrix_is_4 ;matrix_is_5 ;matrix_is_6 ;matrix_is_7 ;;matrix_is_8]);
else
matrix_is=eval( eval_matrix(1:size(eval_matrix,1)*size(eval_matrix,2)));

windowed_data=data([ matrix_is]);
if (step_size-floor(step_size))~=0
      windowed_data=[windowed_data n_data'];
end
end
clear matrix_is*

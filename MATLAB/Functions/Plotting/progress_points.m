function progress_points(n,N)
    if n == 1
        fprintf(append('\nProgress bar: \n|',repmat('-',1,N),'|\n '));
    end
    fprintf('o')
    if n == N
        fprintf('\n')
    end
end
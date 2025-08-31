function output_matrix = convert_adjacency_matrix(adj_matrix)
    [N, ~] = size(adj_matrix);
    output_matrix = zeros(N, N);

    for i = 1:N
        neighbors = find(adj_matrix(i, :));
        
        output_matrix(i, 1:length(neighbors)) = neighbors;
    end
end

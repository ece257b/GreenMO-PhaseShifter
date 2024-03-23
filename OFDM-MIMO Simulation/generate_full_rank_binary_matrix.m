function matrix = generate_full_rank_binary_matrix(prm)
    if(prm.numTxAntenna==prm.numUsers)
        matrix=eye(prm.numTxAntenna);
        return
    end
    while true
        matrix = randi([0 1], prm.numTxAntenna, prm.numUsers);  % Generate random binary matrix
        if rank(matrix) == min(prm.numTxAntenna,prm.numUsers )  % Check if matrix has full rank
            break;
        end
    end
end
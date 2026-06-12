function atom_ave = average_atom_z_direction(atom)
index_up = atom(3,:)<0;
index_down = atom(3,:)>0;
atom_ave = zeros(size(atom));
atom_up = atom(:,index_up);
atom_down = atom(:,index_down);
atom_up_ave = atom_up;
atom_down_ave = atom_down;
for tempi = 1:size(atom_up,2)
    dif=(atom_up-atom_up(:,tempi));
    dis=sqrt(sum(dif.^2,1));
    temp_index = dis<2.6;
    atom_up_ave(3,tempi) = mean(atom_up(3,temp_index));
    clear temp_index
end
clear tempi
for tempi = 1:size(atom_down,2)
    dif=(atom_down-atom_down(:,tempi));
    dis=sqrt(sum(dif.^2,1));
    temp_index = dis<2.6;
    atom_down_ave(3,tempi) = mean(atom_down(3,temp_index));
    clear temp_index
end
clear tempi
atom_ave(:,index_up) = atom_up_ave;
atom_ave(:,index_down) = atom_down_ave;
end
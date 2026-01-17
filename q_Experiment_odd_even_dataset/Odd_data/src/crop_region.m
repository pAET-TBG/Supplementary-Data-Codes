function [atomup,atomup_sim,labelup,atomdown,atomdown_sim,labeldown,atomupgroup,atomdowngroup] = crop_region(atomup_exp,atomdown_exp,atomup_sim,atomdown_sim,labelup,labeldown,atomupgroup,atomdowngroup)

atomup    = atomup_exp;                                      % store the inital position for distance calculate
atomdown  = atomdown_exp;                                    % store the inital position for distance calculate
clear atom

atom      = [atomup,atomdown];
label     = [labelup,labeldown];
atomgroup = [atomupgroup,atomdowngroup];

atomup(1,:) = atomup(1,:) - min(atom(1,:)) + 10;
atomup(2,:) = atomup(2,:) - min(atom(2,:)) + 10;
atomup(3,:) = atomup(3,:) - mean(atom(3,:));

atomup_sim(1,:) = atomup_sim(1,:) - min(atom(1,:)) + 10;
atomup_sim(2,:) = atomup_sim(2,:) - min(atom(2,:)) + 10;
atomup_sim(3,:) = atomup_sim(3,:) - mean(atom(3,:));

atomdown(1,:) = atomdown(1,:) - min(atom(1,:)) + 10;
atomdown(2,:) = atomdown(2,:) - min(atom(2,:)) + 10;
atomdown(3,:) = atomdown(3,:) - mean(atom(3,:));

atomdown_sim(1,:) = atomdown_sim(1,:) - min(atom(1,:)) + 10;
atomdown_sim(2,:) = atomdown_sim(2,:) - min(atom(2,:)) + 10;
atomdown_sim(3,:) = atomdown_sim(3,:) - mean(atom(3,:));

atom(1,:) = atom(1,:) - min(atom(1,:)) + 10;
atom(2,:) = atom(2,:) - min(atom(2,:)) + 10;

%% Sub-region 
points    = [3000 11400
             6000 12220
             7700 9570
             7700 2400 
             6400 130
             3400 130
             1900 2500
             1600 8500]';  

%% crop the subregion
atom      = [atomup,atomdown];
atom_sim  = [atomup_sim,atomdown_sim];
isInside  = inpolygon(atom(1,:), atom(2,:), points(1,:), points(2,:));

atom      = atom(:,isInside);
atom_sim  = atom_sim(:,isInside);
label     = label(isInside);
atomgroup = atomgroup(isInside);

atomup      = atom(:,atom(3,:)<0);
atomup_sim  = atom_sim(:,atom_sim(3,:)<0); 
labelup     = label(atom(3,:)<0); 
atomupgroup = atomgroup(atom(3,:)<0);

atomdown      = atom(:,atom(3,:)>0);
atomdown_sim  = atom_sim(:,atom_sim(3,:)>0); 
labeldown     = label(atom_sim(3,:)>0); 
atomdowngroup = atomgroup(atom(3,:)>0);
clear isInside

end



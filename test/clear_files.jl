
function clean_directory()
    rm_folders = ["Initial_Conditions", "Results", "Results_vtu"]
    rm.(rm_folders;force=true, recursive=true)
end

clean_directory()
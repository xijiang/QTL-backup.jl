"""
Return the size of available memory on a Linux alike system.
Or, the free memory size is returned.
""" 
function memavail()
    if Sys.isunix()
        for line in eachline("/proc/meminfo")
            if occursin("MemAv", line)
                return parse(Int, split(line)[2]) * 1024
            end
        end
    else
        return Sys.free_memory()
    end
end

            

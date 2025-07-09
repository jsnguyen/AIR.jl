function autolog(f, log_filename)

    if splitext(log_filename)[end] != ".log"
        error("Logfile must have a .log extension")
    end

    dirname(log_filename) |> x -> if (!isdir(x) || !isempty(x)) mkpath(x) end

    open(log_filename, "w") do logfile
        console_logger = ConsoleLogger(stdout)
        file_logger = SimpleLogger(logfile)
        demux_logger = TeeLogger(console_logger, file_logger)
        global_logger(demux_logger)
        f()
    end
end

# needs to be a macro to get the file name at compile time
macro logname()
    sf = String(__source__.file)
    filename = basename(sf)
    filedir = dirname(sf)
    log_path = joinpath(filedir, "logs", "$(filename).log")
    return esc(log_path)
end

macro autolog(block)
    sf = String(__source__.file)
    filename = basename(sf)
    filedir = dirname(sf)
    log_path = joinpath(filedir, "logs", "$(filename).log")

    return esc(quote
            autolog($log_path) do
                $block
            end
        end)
end

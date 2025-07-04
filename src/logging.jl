function autolog(f, log_filename)

    if splitext(log_filename)[end] != ".log"
        error("Logfile must have a .log extension")
    end

    open(log_filename, "w") do logfile
        console_logger = ConsoleLogger(stdout)
        file_logger = SimpleLogger(logfile)
        demux_logger = TeeLogger(console_logger, file_logger)
        global_logger(demux_logger)
        f()
    end
end
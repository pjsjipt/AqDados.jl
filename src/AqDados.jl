module AqDados

export AqDadosFile, numchans, numsamples, samplingrate, samplingtimes
export chanunits, channames, daqtime

using Dates


function read_str(io, buflen)
    n = Int(read(io, UInt8))
    buf = read(io, buflen)
    return String(buf[1:n])
end

    

seekrel(io, n) = seek(io, position(io)+n)


struct AqDadosChan
    idx::Int
    name::String
    unit::String
    lsup::Float32
    linf::Float32
    used::Bool
end

struct AqDadosTEMHeader
    date::DateTime
    comment::String
    nc::Int
    fs::Float32
    ns::Int32
    chans::Vector{AqDadosChan}
    ev::Vector{UInt32}
    excomm::Vector{String}
    fileoffs::Int
    sampleformat::DataType
    samplesize::Int
    binrange::Int
end

    
    
function fileheader_tem(io)

    aqtime = read!(io, zeros(UInt16, 7))
    date = DateTime(aqtime[3], aqtime[2], aqtime[1],
                    aqtime[4], aqtime[5], aqtime[6], aqtime[7])
    
    comment = read_str(io, 30)
    # Read channel data

    chans = AqDadosChan[]
    
    ns = 0
    for i in 1:32
        # Ler o nome do canal
        aname = read_str(io, 15)
        aunit = read_str(io, 7)
        asup = read(io, Float32)
        ainf = read(io, Float32)
        used = read(io, UInt8) != 0
        seekrel(io, 2)
        if used
            ns += 1
            push!(chans, AqDadosChan(ns, aname, aunit, asup, ainf, used))
        end
        
    end

    
    # Get the number of aquired channels
    nc = read(io, UInt16)

    
    # Get the sampling frequency
    seekrel(io, 4)
    fs = read(io, Float32)
    ns = read(io, Int32)

    seekrel(io, 252)

    # Read extra comments
    excom = String[]
    for i in 1:6
        push!(excom, read_str(io, 30))
    end

    # Ge the number of events:
    ne = Int(read(io, UInt16))

    # Get the array of events
    ev = read!(io, zeros(UInt32, ne))

    header = AqDadosTEMHeader(date, comment, nc, fs, ns, chans, ev, excom,
                              1747, Int16, 2, 65536)
    return header
end


function aqdados_data(io, header)

    seek(io, header.fileoffs)

    nc = header.nc
    ns = header.ns
    Xi = read!(io, zeros(header.sampleformat, nc, ns))

    lsup = [ch.lsup for ch in header.chans]
    linf = [ch.linf for ch in header.chans]
    m = (lsup .- linf) ./ header.binrange
    b = m .* header.binrange ./ 2 .+ linf
    println(m)
    println(b)
    
    X = zeros(Float64, ns, nc)

    for ich in 1:header.nc
        mc = m[ich]
        bc = b[ich]
        for it in 1:header.ns
            X[it,ich] = mc * Xi[ich,it] + bc
        end
    end
    
    return X
    
    
end

function aqdados_tem_file(fname)
    open(fname, "r") do io
        header = fileheader_tem(io)
        X = aqdados_data(io, header)
        header, X
    end
    
end

struct AqDadosFile{T}
    header::AqDadosTEMHeader
    data::Matrix{T}
end

function AqDadosFile(io::IO)
    seek(io, 0)
    header = fileheader_tem(io)
    data = aqdados_data(io, header)
    AqDadosFile(header, data)
end
AqDadosFile(fname::AbstractString) = open(fname, "r") do io
    AqDadosFile(io)
end

AqDadosFile(bytes::AbstractVector{UInt8}) = AqDadosFile(IOBuffer(bytes))

numchans(d::AqDadosFile) = d.header.nc

channames(d::AqDadosFile) = [ch.name for ch in d.header]
channames(d::AqDadosFile, i) = d.header.chans[i]

chanunits(d::AqDadosFile) = [ch.unit for ch in d.header]
chanunits(d::AqDadosFile, i) = d.header.chans[i].unit

samplingrate(d::AqDadosFile) = d.header.fs
samplingtimes(d::AqDadosFile) = range(0.0, length=d.header.ns,step=1/d.header.fs)

daqtime(d::AqDadosFile) = d.header.date

numsamples(d::AqDadosFile) = d.header.ns

Base.getindex(d::AqDadosFile,i) = d.data[i]
Base.getindex(d::AqDadosFile) = d.data


end

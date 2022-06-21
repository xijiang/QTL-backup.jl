function mytest(sdir, nc, tdir, δ)
    fg   = "$sdir/G.bin"
    fpiv = "$sdir/piv.bin"
    fid  = "$sdir/id.bin"
    
    # fg11, fg21, fd22 = QTL.Mat.extract_subs(fg, fpiv, nc, tdir)
    fg11 = "$tdir/g11-$(nc÷1000)k.bin"
    fg21 = "$tdir/g21-$(nc÷1000)k.bin"
    fd22 = "$tdir/d22-$(nc÷1000)k.bin"
    QTL.Mat.approxgi(fg11, fg21, fd22, fpiv, fid, tdir, δ = δ)
end

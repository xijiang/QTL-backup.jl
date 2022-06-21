function mytest(sdir, nc, tdir, δ)
    fg   = "$sdir/G.bin"
    fpiv = "$sdir/piv.bin"
    fid  = "$sdir/id.bin"
    
    fg11, fg21, fd22 = QTL.Mat.extract_subs(fg, fpiv, nc, tdir)
    QTL.Mat.approxgi(fg11, fg21, fd22, fpiv, fid, tdir, δ = δ)
end

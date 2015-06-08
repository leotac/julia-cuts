function writeULS(path::String, T::Int; seed::Int=3)
   f = open(path,"w")
   srand(seed)
   write(f, "$T\n")
   write(f, "$(reduce(*, map(x -> "$x ", rand(1:4, T))))\n") # c
   write(f, "$(reduce(*, map(x -> "$x ", rand(1:3, T))))\n") # h
   write(f, "$(reduce(*, map(x -> "$x ", rand(5:10, T))))\n") # K
   write(f, "$(reduce(*, map(x -> "$x ", rand(1:10, T))))\n") # d
   close(f)
end

function readULS(path::String)
   f = open(path)
   T = int(readline(f))
   c = map(float,split(strip(readline(f))))
   h = map(float,split(strip(readline(f))))
   K = map(float,split(strip(readline(f))))
   d = map(float,split(strip(readline(f))))
   close(f)
   T, c, h, K, d
end

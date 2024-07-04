using ADerrors
using BDIO
using ALPHAio
import ADerrors: err
using ADjson

a = uwreal(rand(2016), "H101") #uwreal([0.1443, 0.0017], "t0")
b = uwreal(rand(1540), "N300")

c = a * exp(b)
uwerr(c)
details(c)

## saving
fname = tempname()
io = IOBuffer()
fb = ALPHAdobs_create(fname, io)
ALPHAdobs_write(fb, c)
ALPHAdobs_close(fb)

## reading
fb = BDIO_open(fname, "r")
res = 0.0
while ALPHAdobs_next_p(fb)
    d = ALPHAdobs_read_parameters(fb)
    res = ALPHAdobs_read_next(fb)
end
BDIO_close!(fb)

## alternative reading 
fb = BDIO_open(fname, "r")
BDIO_seek!(2)
tt = read_uwreal(fb)
BDIO_close!(fb)
## test derivatives

der_c_a = derivative(c, a) # original variable
der_c_b = derivative(c, b) # original variable

der_res_a = derivative(res, a) # variable read from BDIO
der_res_b = derivative(res, b) # variable read from BDIO


## why?
res.prop == c.prop
res.der == c.der
# these two vectors are not the same in the original variable and the one read from BDIO


## test with json
dump_to_json(c, fname, "test")

check = load_json(fname)
uwerr(check)
check
c
derivative(check,a)
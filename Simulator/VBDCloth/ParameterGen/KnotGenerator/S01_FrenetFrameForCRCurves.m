t = sym('t');
ts = sym('t0', [4,1]);

tangent_norm =  sym('tangent_norm');

p1 = sym('p1', [3,1]);
p2 = sym('p2', [3,1]);
p3 = sym('p3', [3,1]);
p4 = sym('p4', [3,1]);

ps = {p1, p2, p3, p4};


x = cubicCatmullRom(t, ps, ts);



tangent = diff(x,t)

tangent_n = tangent / tangent_norm
normal = diff(tangent_n, t);


normal = simplify(normal)
%fprintf('%s',char(normal));


binormal = cross(tangent_n, normal);


function x = interpolate(x1, x2, t, t0, t1)

     x=x1 * (t1-t)/(t1-t0) + x2 * (t-t0)/(t1-t0);
end

function C12 = cubicCatmullRom(t, ps, ts)
    L01 = interpolate(ps(1), ps(2), t, ts(1), ts(2));
    L12 = interpolate(ps(2), ps(3), t, ts(2), ts(3));
    L23 = interpolate(ps(3), ps(4), t, ts(3), ts(4));

    L012 = interpolate(L01, L12, t, ts(1), ts(3));
    L123 = interpolate(L12, L23, t, ts(2), ts(4));

    C12 = interpolate(L012, L123, t, ts(2), ts(3));

end




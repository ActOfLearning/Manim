from manimlib.imports import *
import numpy as np

class Angle(VGroup):

    CONFIG = {
        'radius': 0.5,
        'color': RED,
        'opacity': 0.4,
        'stroke_width': 5,
        'below_180': True,
    }

    def __init__(self, A, O, B, **kwargs):

        VMobject.__init__(self, **kwargs)
        OA, OB = A-O, B-O
        if self.below_180:
            theta = np.angle(complex(*OA[:2])/complex(*OB[:2])) # angle of OB to OA
        else:
            theta = TAU + np.angle(complex(*OA[:2])/complex(*OB[:2]))

        if abs(abs(theta) - PI / 2) <= 1e-2:
            self.add(Elbow(angle = np.angle(complex(*OB[:2])), color = self.color, width = self.radius).shift(O))
        else:
            self.add(Sector(inner_radius = 0, outer_radius = self.radius, angle = theta, arc_center = O, color = self.color, fill_opacity = self.opacity, start_angle = Line(O, B).get_angle()))
            self.add(Arc(start_angle=Line(O, B).get_angle(), angle=theta, radius=self.radius,
                     stroke_width=self.stroke_width, color=self.color, arc_center=O))

def Circle_Inversion(obj, about_point = ORIGIN, inversion_radius = 1.0):
    def inverse_of_a_point(pt0):
        norm_of_pt = np.dot(pt0, pt0)
        if norm_of_pt == 0:
            norm_of_pt += 0.0001
        return pt0 / norm_of_pt
    def inverse_of_a_circle(obj):
        center_of_the_circle = obj.get_center()
        radius_of_the_circle = np.linalg.norm(center_of_the_circle - obj.get_start())
        inversion_factor = np.dot(center_of_the_circle, center_of_the_circle) - radius_of_the_circle ** 2
        if inversion_factor == 0:
            res_object = Line(10 * UP, 10 * DOWN)
            res_object.move_to(RIGHT * (1 / 2 / radius_of_the_circle))
            rot_angle = np.angle(center_of_the_circle[0] + center_of_the_circle[1] * 1j)
            res_object.rotate_about_origin(rot_angle)
            return res_object
        inversion_factor = 1.0 / inversion_factor
        return Circle(radius = abs(inversion_factor) * radius_of_the_circle).move_to(inversion_factor * center_of_the_circle).set_color(TEAL)
    def circum_center(pt0):
        pt = pt0[:]
        side_a, side_b, side_c = pt[2] - pt[1], pt[0] - pt[2], pt[1] - pt[0]
        side_a, side_b, side_c = np.dot(side_a, side_a), np.dot(side_b, side_b), np.dot(side_c, side_c)
        circum_x, circum_y, circum_z = side_a * (side_b + side_c - side_a), side_b * (side_c + side_a - side_b), side_c * (side_a + side_b - side_c)
        return (circum_x * pt[0] + circum_y * pt[1] + circum_z * pt[2]) / (circum_x + circum_y + circum_z)
    def inverted_arc(c0, r0, ang0):
        return Arc(ang0[0], ang0[1] - ang0[0] + (2 * PI if ang0[1] < ang0[0] else 0), radius = r0, arc_center = c0)
    def arg_of_a_point(pt0):
        return np.angle(pt0[0] + 1j * pt0[1])
    def inverse_of_a_line(obj):
        line_st, line_en = obj.get_start(), obj.get_end()
        inv_line_start = inverse_of_a_point(line_st)
        inv_line_end = inverse_of_a_point(line_en)
        if abs(np.linalg.det([line_st, line_en, OUT])) >= 0.001:
            circum = circum_center([ORIGIN, inv_line_start, inv_line_end])
            circum_norm = np.linalg.norm(circum)
            temp_mat = np.array([[0, -1, 0], [1, 0, 0], [0, 0, 1]])
            temp_mat = np.dot(np.dot(temp_mat, np.transpose(line_st)), line_en)
            return inverted_arc(circum, circum_norm, [arg_of_a_point(inv_line_start - circum), arg_of_a_point(inv_line_end - circum)] if temp_mat > 0 else [arg_of_a_point(inv_line_end - circum), arg_of_a_point(inv_line_start - circum)])
        else:
            if np.dot(line_st, line_en) > 0:
                return Line(inv_line_start, inv_line_end)
            else:
                return VGroup(Line(inv_line_start, 6 * inv_line_start / (0.1 + np.linalg.norm(inv_line_start))), Line(inv_line_end, 6 * inv_line_end / (0.1 + np.linalg.norm(inv_line_end))))
    def inverse_of_a_polygon(obj):
        polygon_vertices = obj.get_points()
        inverted_points = VGroup()
        for i in range(1, len(polygon_vertices)):
            inverted_points.add(inverse_of_a_line(Line(polygon_vertices[i], polygon_vertices[i - 1])))
        return inverted_points.set_color(TEAL)
    shifted_obj = obj.copy()
    shifted_obj.shift(-about_point)
    shifted_obj.scale(1 / inversion_radius, about_point = ORIGIN)
    if type(obj).__name__ == "Circle":
        res_obj = inverse_of_a_circle(shifted_obj)
    elif type(obj).__name__ == "Dot":
        res_obj = Dot(radius = DEFAULT_DOT_RADIUS / inversion_radius).move_to(inverse_of_a_point(shifted_obj.get_center()))
    elif type(obj).__name__ == "VGroup":
        res_obj = VGroup()
        for objl in shifted_obj:
            res_obj.add(Circle_Inversion(objl))
    else:
        res_obj = inverse_of_a_polygon(shifted_obj)
    res_obj.scale(inversion_radius, about_point = ORIGIN)
    res_obj.shift(about_point)
    return res_obj.set_color(TEAL)

def ParallalLine(pt0, line0):
    res = VGroup()
    if type(pt0).__name__ == "Dot" and type(line0).__name__ == "Line":
        line_start, line_end = line0.get_start_and_end()
        p_cent = pt0.get_center()
        if np.linalg.norm(line_end - line_start) > 0:
            normalized_line = (line_end - line_start) / np.linalg.norm(line_end - line_start)
            return Line(p_cent, p_cent + normalized_line)
    return res

def PerpendicularFoot(pt0, line0):
    res = VGroup()
    if type(pt0).__name__ == "Dot" and type(line0).__name__ == "Line":
        line_start, line_end = line0.get_start_and_end()
        p_cent = pt0.get_center()
        if np.linalg.norm(line_end - line_start) > 0:
            normalized_line = (line_end - line_start) / np.linalg.norm(line_end - line_start)
            line_start += np.dot(p_cent - line_start, normalized_line) * normalized_line
            return Dot().move_to(line_start)
    return res

def PerpendicularLine(pt0, line0):
    line_start, line_end = line0.get_start_and_end()
    line_unit = (line_end - line_start) / np.linalg.norm(line_end - line_start)
    perp_unit = np.cross(line_unit, OUT)
    perp_foot = PerpendicularFoot(pt0, line0).get_center()
    return Line(perp_foot, perp_unit + perp_foot)

def DottoCircle(obj1):
    cen1, rad1 = ORIGIN, 0.0
    if type(obj1).__name__ == "Dot" or type(obj1).__name__ == "Circle":
        cen1 = obj1.get_center()
        if type(obj1).__name__ == "Dot":
            rad1 = 0.0
        else:
            rad1 = np.linalg.norm(cen1 - obj1.get_start())
    return cen1, rad1

def RadicalAxis(obj1, obj2):
    if type(obj1).__name__ == "Dot" or type(obj1).__name__ == "Circle":
        cen1, rad1 = DottoCircle(obj1)
    else:
        print(type(obj1).__name__)
    if type(obj2).__name__ == "Dot" or type(obj2).__name__ == "Circle":
        cen2, rad2 = DottoCircle(obj2)
    else:
        print(type(obj2).__name__)
    dist_cent = np.linalg.norm(cen1 - cen2)
    if dist_cent > 0:
        res = rad1 ** 2 - rad2 ** 2
        res = dist_cent + res / dist_cent
        res /= 2
        res /= dist_cent
        res = cen1 + res * (cen2 - cen1)
        return PerpendicularLine(Dot().move_to(res), Line(cen1, cen2))

def LineSegIntersection(obj1, obj2):
    res_obj = VGroup()
    obj1_start, obj1_end = obj1.get_start_and_end()
    a = obj1_end - obj1_start
    obj2_start, obj2_end = obj2.get_start_and_end()
    b = obj2_end - obj2_start
    c = obj2_start - obj1_start
    acrossb = np.cross(a, b)
    bcrossa = np.cross(b, a)
    acrossc = np.cross(a, c)
    ccrossb = np.cross(c, b)
    acrossb_dot = np.dot(acrossb, acrossb)
    bcrossa_dot = np.dot(bcrossa, bcrossa)
    if acrossb_dot > 0:
        s = np.dot(ccrossb, acrossb) / acrossb_dot
        t = np.dot(acrossc, bcrossa) / bcrossa_dot
        if 0 <= s <= 1 and 0 <= t <= 1:
            res_obj.add(Dot(obj1_start + s * a))
    return res_obj

def RayIntersection(obj1, obj2):
    #needs more work
    res_obj = VGroup()
    obj1_start, obj1_end = obj1.get_start_and_end()
    ab = obj1_end - obj1_start
    ab_unit = ab / np.linalg.norm(ab)
    obj2_start, obj2_end = obj2.get_start_and_end()
    cd = obj2_end - obj2_start
    cd_unit = cd / np.linalg.norm(cd)
    ac = obj2_start - obj1_start
    v = np.cross(ab_unit, cd_unit)
    if np.linalg.norm(v) > 0:
        v_unit = v / np.linalg.norm(v)
        v_sq = np.dot(v, v)
        s1 = np.dot(ac, np.cross(cd_unit, v_unit)) / v_sq
        s2 = np.dot(ac, np.cross(ab_unit, v_unit)) / v_sq
        print(s1, s2, obj2_start, obj2_end)
        if 0 < s1 < 1 and 0 < s2 < 1:
            return (obj1_start + ab_unit * s1 + obj2_start + cd_unit * s2) / 2
    return res_obj

def Intersection(*args):
    res_obj = VGroup()
    obj1, obj2 = args[0], args[1]
    if type(obj1).__name__ == "Line":
        obj1_start, obj1_end = obj1.get_start_and_end()
        ab = obj1_end - obj1_start
        ab_unit = ab / np.linalg.norm(ab)
        if type(obj2).__name__ == "Line":
            obj2_start, obj2_end = obj2.get_start_and_end()
            cd = obj2_end - obj2_start
            ac = obj2_start - obj1_start
            abcd = np.cross(ab, cd)
            abcd_n = np.dot(abcd, abcd)
            if abcd_n != 0:
                abcd = abcd / abcd_n
                abcd = np.dot(np.cross(ac, cd), abcd)
                abcd = obj1_start + ab * abcd
                res_obj.add(Dot(abcd))
        elif type(obj2).__name__ == "Circle":
            cen = Dot().move_to(obj2.get_center())
            p_foot = PerpendicularFoot(cen, obj1)
            p_dist = np.linalg.norm(cen.get_center() - p_foot.get_center())
            rad = np.linalg.norm(obj2.get_center() - obj2.get_start())
            if p_dist == rad:
                res_obj.add(p_foot)
            elif p_dist < rad:
                chord_dist = rad * rad - p_dist * p_dist
                chord_dist = np.sqrt(chord_dist)
                res_obj.add(Dot().move_to(p_foot.get_center() + chord_dist * ab_unit), Dot().move_to(p_foot.get_center() - chord_dist * ab_unit))
    elif type(obj1).__name__ == "Circle":
        cen1, cen2 = obj1.get_center(), obj2.get_center()
        dist_cen = np.linalg.norm(cen1 - cen2)
        rad1, rad2 = np.linalg.norm(cen1 - obj1.get_start()), np.linalg.norm(cen2 - obj2.get_start())
        if abs(dist_cen - rad1 - rad2) < 1e-3:
            res_obj.add(Dot().move_to((rad2 * cen1 + rad1 * cen2) / (rad1 + rad2)))
        elif abs(dist_cen - rad1 + rad2) < 1e-3:
            res_obj.add(Dot().move_to((-rad2 * cen1 + rad1 * cen2) / (rad1 - rad2)))
        elif abs(dist_cen + rad1 - rad2) < 1e-3:
            res_obj.add(Dot().move_to((rad2 * cen1 - rad1 * cen2) / (-rad1 + rad2)))
        elif dist_cen < rad1 + rad2:
            p_foot = RadicalAxis(obj1, obj2)
            #p_foot = PerpendicularLine(p_foot, Line(cen1, cen2))
            return Intersection(p_foot, obj1)
    return res_obj

def Midpoint(*args):
    if type(args[0]).__name__ == "Line":
        line_start, line_end = args[0].get_start_and_end()
        mid_pt = (line_start + line_end) / 2
        mid_pt = Dot().move_to(mid_pt)
        return mid_pt
    elif type(args[0]).__name__ == "Dot" and type(args[1]).__name__ == "Dot":
        mid_pt = (args[0].get_center() + args[1].get_center()) / 2
        mid_pt = Dot().move_to(mid_pt)
        return mid_pt

def PerpendicularBisector(*args):
    mid_pt = Midpoint(*args)
    if type(args[0]).__name__ == "Line":
        print("It's a line..")
        return PerpendicularLine(mid_pt, args[0])
    elif type(args[0]).__name__ == "Dot" and type(args[1]).__name__ == "Dot":
        connec_line = Line(args[0].get_center(), args[1].get_center())
        return PerpendicularLine(mid_pt, connec_line)

def AngularBisector(*args):
    if type(args[0]).__name__ == "Line" and type(args[1]).__name__ == "Line":
        intersec = Intersection(*args)
        if len(list(intersec)) > 0:
            return AngularBisector(Dot().move_to(args[0].get_end()), intersec, Dot().move_to(args[1].get_end()))
        else:
            return VGroup()
    elif type(args[0]).__name__ == "Dot" and type(args[1]).__name__ == "Dot" and type(args[2]).__name__ == "Dot":
        pt_a = args[0].get_center()
        pt_b = args[1].get_center()
        pt_c = args[2].get_center()
        vec_ba = pt_a - pt_b
        vec_bc = pt_c - pt_b
        return Line(pt_b, vec_ba + vec_bc)

def RadicalCenter(*args):
    if type(args[0]).__name__ == "Circle" and type(args[1]).__name__ == "Circle" and type(args[2]).__name__ == "Circle":
        radic12 = RadicalAxis(args[0], args[1])
        radic23 = RadicalAxis(args[1], args[2])
        return Intersection(radic12, radic23)
    return Dot()

def HomotheticCenter(obj1, obj2):
    cen1, cen2 = obj1.get_center(), obj2.get_center()
    rad1, rad2 = np.linalg.norm(cen1 - obj1.get_start()), np.linalg.norm(cen2 - obj2.get_start())
    res_obj = VGroup()
    res_obj.add(Dot().move_to((rad1 * cen2 + rad2 * cen1) / (rad1 + rad2)))
    if rad1 != rad2:
        res_obj.add(Dot().move_to((rad1 * cen2 - rad2 * cen1) / (rad1 - rad2)))
    return res_obj

def CircleTangentLines(obj, pt):
    cen = obj.get_center()
    rad = np.linalg.norm(cen - obj.get_start())
    dist = np.linalg.norm(cen - pt.get_center())
    res_obj = VGroup()
    if dist < rad:
        temp_c = Intersection(PerpendicularLine(pt, Line(obj.get_center(), pt.get_center())), obj)
        for i in temp_c:
            res_obj.add(Line(pt.get_center(), i.get_center()))
    elif dist == rad:
        res_obj.add(Dot(pt.get_center()))
    else:
        temp_c = Circle(arc_center = (cen + pt.get_center()) / 2, radius = dist / 2)
        temp_c = Intersection(temp_c, obj)
        for i in temp_c:
            res_obj.add(Line(pt.get_center(), i.get_center()))
    return res_obj

def ThreePointCircle(ptA, ptB, ptC):
    cen = Intersection(
        PerpendicularBisector(ptA, ptB),
        PerpendicularBisector(ptB, ptC)
    )
    if len(cen) <= 0 or Line(ORIGIN, cen[0].get_center()).get_length() > 1000:
        temp = PerpendicularLine(ptB, Line(ptA.get_center(), ptC.get_center()))
        temp = Line(temp.get_start(), temp.get_start() + 10 ** 5 * (temp.get_end() - temp.get_start()))
        return Circle(arc_center = temp.get_end(), radius = temp.get_length())
    cen = cen[0]
    return Circle(arc_center = cen.get_center(), radius = Line(cen.get_center(), ptA.get_center()).get_length())

def ToCircle(obj):
    if type(obj).__name__ == "Circle":
        return obj
    elif type(obj).__name__ == "Dot":
        return Circle(arc_center = obj.get_center(), radius = 0.001)
    elif type(obj).__name__ == "Line":
        temp = PerpendicularLine(Dot(), obj)
        temp = Line(temp.get_start(), temp.get_start() + 10 ** 5 * (temp.get_end() - temp.get_start()))
        return Circle(arc_center = temp.get_end(), radius = temp.get_length())
    else:
        print("You Idiot.. Dont goof around... Gimme a Dot, Line or a Circle...")
        return Circle()

def ApolloniusCircles(obj1, obj2, obj3):
    ca, cb, cc = ToCircle(obj1), ToCircle(obj2), ToCircle(obj3)
    radcen = RadicalCenter(ca, cb, cc)
    abhcen = HomotheticCenter(ca, cb)
    bchcen = HomotheticCenter(cb, cc)
    cahcen = HomotheticCenter(cc, ca)
    haxis = VGroup(
        Line(abhcen[0].get_center(), bchcen[0].get_center()),
        Line(bchcen[0].get_center(), cahcen[0].get_center()),
        Line(cahcen[0].get_center(), abhcen[0].get_center())
    )
    if len(abhcen) > 1 and len(bchcen) > 1:
        haxis.add(Line(abhcen[-1].get_center(), bchcen[-1].get_center()))
    elif len(bchcen) > 1 and len(cahcen) > 1:
        haxis.add(Line(bchcen[-1].get_center(), cahcen[-1].get_center()))
    elif len(cahcen) > 1 and len(abhcen) > 1:
        haxis.add(Line(cahcen[-1].get_center(), abhcen[-1].get_center()))
    solns = VGroup()
    for obj in haxis:
        tpts = VGroup()
        for circ in [ca, cb, cc]:
            irad = Line(circ.get_center(), circ.get_start()).get_length()
            perpft = PerpendicularFoot(Dot(circ.get_center()), obj)
            perpft = Circle_Inversion(perpft, about_point = circ.get_center(), inversion_radius = irad)
            perpft = Line(radcen.get_center(), perpft.get_center())
            perpft = Intersection(perpft, circ)
            tpts.add(perpft)
        if len(tpts[0]) < 2 or len(tpts[1]) < 2 or len(tpts[2]) < 2:
            continue
        for j in range(4):
            tcirc = ThreePointCircle(tpts[0][0], tpts[1][j % 2], tpts[2][(j // 2) % 2])
            tcirccen = tcirc.get_center()
            tcircrad = Line(tcirccen, tcirc.get_start()).get_length()
            flag = 0
            for ccir in [ca, cb, cc]:
                cencenlen = Line(tcirccen, ccir.get_center()).get_length()
                ccirrad = Line(ccir.get_center(), ccir.get_start()).get_length()
                if abs(cencenlen - (ccirrad + tcircrad)) <= 0.001 or abs(cencenlen - (ccirrad - tcircrad)) <= 0.001 or abs(cencenlen - (-ccirrad + tcircrad)) <= 0.001:
                    continue
                else:
                    flag = 1
                    break
            if flag == 0:
                solns.add(tcirc)
                tcirc = ThreePointCircle(tpts[0][1], tpts[1][1 - (j % 2)], tpts[2][1 - (j // 2) % 2])
                solns.add(tcirc)
                break
    solns.set_color(GREEN)
    return solns

class IntroTest(GraphScene, MovingCameraScene):
    CONFIG = {
        "x_min": -10, "x_max": 10, "x_axis_width": 20,
        "y_min": -5, "y_max": 5, "y_axis_height": 10,
        "axes_color": BLACK,
        "graph_origin": ORIGIN,
    }
    def setup(self):
        GraphScene.setup(self)
        MovingCameraScene.setup(self)
    
    def construct(self):
        self.setup_axes()
        self.play(FadeOut(self.axes))
        self.wait()
        nplane = VGroup(DoubleArrow(8 * LEFT, 8 * RIGHT), DoubleArrow(5.5 * DOWN, 5.5 * UP)).set_color(GREY)
        self.play(Write(nplane))
        self.wait()
        planegrid = VGroup()
        for i in range(-10, 10 + 1):
            planegrid.add(Line(5 * UP + i * RIGHT, 5 * DOWN + i * RIGHT))
        for i in range(-5, 5 + 1):
            planegrid.add(Line(i * UP + 10 * RIGHT, i * UP + 10 * LEFT))
        planegrid.set_color(GREY)
        planegrid.fade(0.875)
        self.play(ShowCreation(planegrid), run_time = 3)
        self.wait()
        xa, ya, ra = ValueTracker(2), ValueTracker(2), ValueTracker(1)
        xb, yb, rb = ValueTracker(-2), ValueTracker(2), ValueTracker(1)
        xc, yc, rc = ValueTracker(0), ValueTracker(-2), ValueTracker(0.999)
        def circcreator(aco, bco, cco):
            x1, y1, r1 = aco[0], aco[1], aco[2]
            x2, y2, r2 = bco[0], bco[1], bco[2]
            x3, y3, r3 = cco[0], cco[1], cco[2]
            ca = Circle(arc_center = x1 * RIGHT + y1 * UP, radius = r1)
            cb = Circle(arc_center = x2 * RIGHT + y2 * UP, radius = r2)
            cc = Circle(arc_center = x3 * RIGHT + y3 * UP, radius = r3)
            solns = ApolloniusCircles(ca, cb, cc)
            return VGroup(VGroup(ca, cb, cc), solns)
        def circupdater(obj):
            aco = [xa.get_value(), ya.get_value(), ra.get_value()]
            bco = [xb.get_value(), yb.get_value(), rb.get_value()]
            cco = [xc.get_value(), yc.get_value(), rc.get_value()]
            obj.become(circcreator(aco, bco, cco))
        solgrp = circcreator([xa.get_value(), ya.get_value(), ra.get_value()], [xb.get_value(), yb.get_value(), rb.get_value()], [xc.get_value(), yc.get_value(), rc.get_value()])
        self.play(Write(solgrp[0]))
        for obj in solgrp[1]:
            self.play(Write(obj))
        self.play(
            #xa.increment_value, -1.5,
            #xb.increment_value, -1,
            #xc.increment_value, -2,
            ya.increment_value, -4,
            yb.increment_value, -4,
            yc.increment_value, 4,
            ra.increment_value, 1,
            rb.increment_value, 1,
            rc.increment_value, 1,
            UpdateFromFunc(solgrp, circupdater),
            run_time = 5)
        self.wait(5)

class Test(GraphScene, MovingCameraScene):
    CONFIG = {
        "x_min": -10, "x_max": 10, "x_axis_width": 20,
        "y_min": -5, "y_max": 5, "y_axis_height": 10,
        "axes_color": BLACK,
        "graph_origin": ORIGIN,
    }
    def setup(self):
        GraphScene.setup(self)
        MovingCameraScene.setup(self)
    
    def construct(self):
        self.setup_axes()
        self.play(FadeOut(self.axes))
        self.wait()
        nplane = VGroup(DoubleArrow(8 * LEFT, 8 * RIGHT), DoubleArrow(5.5 * DOWN, 5.5 * UP)).set_color(GREY)
        self.play(Write(nplane))
        self.wait()
        planegrid = VGroup()
        for i in range(-10, 10 + 1):
            planegrid.add(Line(5 * UP + i * RIGHT, 5 * DOWN + i * RIGHT))
        for i in range(-5, 5 + 1):
            planegrid.add(Line(i * UP + 10 * RIGHT, i * UP + 10 * LEFT))
        planegrid.set_color(GREY)
        planegrid.fade(0.875)
        self.play(ShowCreation(planegrid), run_time = 3)
        self.wait()
        '''xa, ya, ra = 0, 0, 1
        xb, yb, rb = 1, 0, 1
        xc, yc, rc = 0.5, 3 ** 0.5 / 2, 1.001'''
        xa, ya, ra = 0, 0, 3.5
        xb, yb, rb = 1, 0, 0.5
        xc, yc, rc = -1, -1, 0.75
        ca = Circle(arc_center = xa * RIGHT + ya * UP, radius = ra)
        cb = Circle(arc_center = xb * RIGHT + yb * UP, radius = rb)
        cc = Circle(arc_center = xc * RIGHT + yc * UP, radius = rc)
        self.play(
            Write(ca), Write(cb), Write(cc)
        )
        radcen = RadicalCenter(ca, cb, cc)
        abhcen = HomotheticCenter(ca, cb)
        bchcen = HomotheticCenter(cb, cc)
        cahcen = HomotheticCenter(cc, ca)
        haxis = VGroup(
            Line(abhcen[0].get_center(), bchcen[0].get_center()),
            Line(bchcen[0].get_center(), cahcen[0].get_center()),
            Line(cahcen[0].get_center(), abhcen[0].get_center())
        )
        if len(abhcen) > 1 and len(bchcen) > 1:
            haxis.add(Line(abhcen[-1].get_center(), bchcen[-1].get_center()))
        elif len(bchcen) > 1 and len(cahcen) > 1:
            haxis.add(Line(bchcen[-1].get_center(), cahcen[-1].get_center()))
        elif len(cahcen) > 1 and len(abhcen) > 1:
            haxis.add(Line(cahcen[-1].get_center(), abhcen[-1].get_center()))
        solns = VGroup()
        #self.play(Write(radcen))
        for obj in haxis:
            tpts = VGroup()
            #self.play(Write(obj))
            for circ in [ca, cb, cc]:
                irad = Line(circ.get_center(), circ.get_start()).get_length()
                perpft = PerpendicularFoot(Dot(circ.get_center()), obj)
                #self.play(Write(perpft))
                perpft = Circle_Inversion(perpft, about_point = circ.get_center(), inversion_radius = irad)
                #self.play(Write(perpft))
                perpft = Line(radcen.get_center(), perpft.get_center())
                #self.play(Write(perpft))
                perpft = Intersection(perpft, circ)
                #self.play(Write(perpft))
                tpts.add(perpft)
            if len(tpts[0]) < 2 or len(tpts[1]) < 2 or len(tpts[2]) < 2:
                continue
            for j in range(4):
                tcirc = ThreePointCircle(tpts[0][0], tpts[1][j % 2], tpts[2][(j // 2) % 2])
                #self.play(Write(tcirc))
                tcirccen = tcirc.get_center()
                tcircrad = Line(tcirccen, tcirc.get_start()).get_length()
                flag = 0
                for ccir in [ca, cb, cc]:
                    cencenlen = Line(tcirccen, ccir.get_center()).get_length()
                    ccirrad = Line(ccir.get_center(), ccir.get_start()).get_length()
                    if abs(cencenlen - (ccirrad + tcircrad)) <= 0.001 or abs(cencenlen - (ccirrad - tcircrad)) <= 0.001 or abs(cencenlen - (-ccirrad + tcircrad)) <= 0.001:
                        continue
                    else:
                        flag = 1
                        break
                if flag == 0:
                    solns.add(tcirc)
                    tcirc = ThreePointCircle(tpts[0][1], tpts[1][1 - (j % 2)], tpts[2][1 - (j // 2) % 2])
                    solns.add(tcirc)
                    break
        solns.set_color(GREEN)
        for obj in solns:
            self.play(Write(obj))
        self.wait(5)

class Intro(GraphScene, MovingCameraScene):
    CONFIG = {
        "x_min": -10, "x_max": 10, "x_axis_width": 20,
        "y_min": -5, "y_max": 5, "y_axis_height": 10,
        "axes_color": BLACK,
        "graph_origin": ORIGIN,
    }
    def setup(self):
        GraphScene.setup(self)
        MovingCameraScene.setup(self)
    
    def construct(self):
        self.setup_axes()
        self.play(FadeOut(self.axes))
        #self.wait()
        #nplane = VGroup(DoubleArrow(8 * LEFT, 8 * RIGHT), DoubleArrow(5.5 * DOWN, 5.5 * UP)).set_color(GREY)
        #self.play(Write(nplane))
        self.wait()
        planegrid = VGroup()
        for i in range(-10, 10 + 1):
            planegrid.add(Line(5 * UP + i * RIGHT, 5 * DOWN + i * RIGHT))
        for i in range(-5, 5 + 1):
            planegrid.add(Line(i * UP + 10 * RIGHT, i * UP + 10 * LEFT))
        planegrid.set_color(GREY)
        planegrid.fade(0.875)
        self.play(ShowCreation(planegrid), run_time = 3)
        self.wait()
        a1, a2, a3 = ValueTracker(0), ValueTracker(0), ValueTracker(0)
        r1, r2, r3 = ValueTracker(1), ValueTracker(0.99), ValueTracker(1)
        p1 = Arc(radius = 2, angle = 120 * DEGREES)
        p2 = Arc(radius = 2, start_angle = 120 * DEGREES, angle = 240 * DEGREES)
        p3 = Arc(radius = 2, start_angle = 240 * DEGREES, angle = 360 * DEGREES)
        c1 = Circle(arc_center = p1.point_from_proportion(a1.get_value()), radius = r1.get_value())
        c2 = Circle(arc_center = p2.point_from_proportion(a2.get_value()), radius = r2.get_value())
        c3 = Circle(arc_center = p3.point_from_proportion(a3.get_value()), radius = r3.get_value())
        soln = ApolloniusCircles(c1, c2, c3)
        #self.play(Write(p1), Write(p2), Write(p3))
        self.play(Write(c1), Write(c2), Write(c3))
        for obj in soln:
            self.add(obj)
            self.wait(0.5)
            self.remove(obj)
        self.wait()
        soln = VGroup(soln[1])
        self.play(Write(soln))
        dyngrp = VGroup(c1, c2, c3, soln)
        def solupd(obj):
            #tc1, tc2, tc3, ts = obj
            cc1 = Circle(arc_center = p1.point_from_proportion(a1.get_value()), radius = r1.get_value())
            cc2 = Circle(arc_center = p2.point_from_proportion(a2.get_value()), radius = r2.get_value())
            cc3 = Circle(arc_center = p3.point_from_proportion(a3.get_value()), radius = r3.get_value())
            temp = ApolloniusCircles(cc1, cc2, cc3)
            temp = VGroup(temp[1]) if len(temp) > 1 else VGroup()
            obj.become(VGroup(cc1, cc2, cc3, temp))
        dyngrp.add_updater(solupd)
        self.add(dyngrp)
        self.play(
            r1.increment_value, -0.975,
            run_time = 3
        )
        self.wait()
        dyngrp.clear_updaters()
        x1 =  ValueTracker(2)
        def solupd1(obj):
            #tc1, tc2, tc3, ts = obj
            cc1 = Circle(arc_center = x1.get_value() * RIGHT, radius = r1.get_value())
            cc2 = Circle(arc_center = p2.point_from_proportion(a2.get_value()), radius = r2.get_value())
            cc3 = Circle(arc_center = p3.point_from_proportion(a3.get_value()), radius = r3.get_value())
            temp = ApolloniusCircles(cc1, cc2, cc3)
            temp = VGroup(temp[1]) if len(temp) > 1 else VGroup()
            obj.become(VGroup(cc1, cc2, cc3, temp))
        dyngrp.add_updater(solupd1)
        self.add(dyngrp)
        self.play(
            r1.increment_value, 150,
            x1.increment_value, 150,
            run_time = 3
        )
        self.wait()
        self.play(
            r1.set_value, 1,
            x1.set_value, p1.point_from_proportion(a1.get_value())[0],
            run_time = 3
        )
        dyngrp.clear_updaters()
        self.wait()
        self.play(FadeOut(soln))
        soln = ApolloniusCircles(c1, c2, c3)
        self.play(FadeIn(soln))
        self.wait()
        dyngrp = VGroup(c1, c2, c3, soln)
        def solupd2(obj):
            #tc1, tc2, tc3, ts = obj
            cc1 = Circle(arc_center = p1.point_from_proportion(a1.get_value()), radius = r1.get_value())
            cc2 = Circle(arc_center = p2.point_from_proportion(a2.get_value()), radius = r2.get_value())
            cc3 = Circle(arc_center = p3.point_from_proportion(a3.get_value()), radius = r3.get_value())
            temp = ApolloniusCircles(cc1, cc2, cc3)
            #temp = VGroup(temp[1]) if len(temp) > 1 else VGroup()
            obj.become(VGroup(cc1, cc2, cc3, temp))
        dyngrp.add_updater(solupd2)
        self.add(dyngrp)
        self.play(
            a1.set_value, 0.333,
            a2.set_value, 0.5,
            a3.set_value, 0.333,
            r1.increment_value, 0.25,
            r3.increment_value, -0.25,
            UpdateFromFunc(dyngrp, solupd2),
            run_time = 3
        )
        self.play(
            a1.set_value, 0.666,
            a2.set_value, 0.666,
            a3.set_value, 0.666,
            r1.increment_value, 0.25,
            r3.increment_value, -0.25,
            UpdateFromFunc(dyngrp, solupd2),
            run_time = 3
        )
        self.play(
            a1.set_value, 0,
            a2.set_value, 0,
            a3.set_value, 0,
            r1.set_value, 1,
            r3.set_value, 1,
            UpdateFromFunc(dyngrp, solupd2),
            run_time = 3
        )
        self.wait()
        dyngrp.clear_updaters()
        gerg = ImageMobject("Gergonne.jpeg").scale(2.5)
        gerg.shift(UP)
        gername = VGroup(
            TextMobject("Joseph Diez Gergonne"),
            TextMobject("(1771 - 1859)"),
        )
        gername.arrange(DOWN)
        gername.add_background_rectangle()
        gername.next_to(gerg, DOWN)
        self.play(
            FadeIn(gerg),
            Write(gername)
        )
        self.wait(5)

class KeyIdeas(GraphScene, MovingCameraScene):
    CONFIG = {
        "x_min": -10, "x_max": 10, "x_axis_width": 20,
        "y_min": -5, "y_max": 5, "y_axis_height": 10,
        "axes_color": BLACK,
        "graph_origin": ORIGIN,
    }
    def setup(self):
        GraphScene.setup(self)
        MovingCameraScene.setup(self)
    
    def construct(self):
        self.setup_axes()
        self.play(FadeOut(self.axes))
        #self.wait()
        #nplane = VGroup(DoubleArrow(8 * LEFT, 8 * RIGHT), DoubleArrow(5.5 * DOWN, 5.5 * UP)).set_color(GREY)
        #self.play(Write(nplane))
        self.wait()
        planegrid = VGroup()
        for i in range(-10, 10 + 1):
            planegrid.add(Line(5 * UP + i * RIGHT, 5 * DOWN + i * RIGHT))
        for i in range(-5, 5 + 1):
            planegrid.add(Line(i * UP + 10 * RIGHT, i * UP + 10 * LEFT))
        planegrid.set_color(GREY)
        planegrid.fade(0.875)
        self.add(planegrid)
        self.wait()
        keyideas = TextMobject("$\\underline{\\text{Key Ideas}}$").scale(2)
        self.play(Write(keyideas))
        self.play(keyideas.shift, 3 * UP)
        self.wait()
        idea1 = TextMobject("$\\bullet$ Power of a point theorem").scale(1.25)
        idea1[0][1:14].set_color(BLUE)
        idea1.move_to(6.125 * LEFT + UP + idea1.get_center() - idea1.get_critical_point(LEFT))
        idea2 = TextMobject("$\\bullet$ Reciprocal relation between poles and polars").scale(1.25)
        idea2[0][1:11].set_color(PINK)
        idea2[0][26:31].set_color(YELLOW)
        idea2[0][34:].set_color(YELLOW)
        idea2.move_to(6.125 * LEFT + DOWN + idea2.get_center() - idea2.get_critical_point(LEFT))
        self.play(Write(idea1), Write(idea2))
        self.wait(5)

class Powerofpoint(GraphScene, MovingCameraScene):
    CONFIG = {
        "x_min": -10, "x_max": 10, "x_axis_width": 20,
        "y_min": -5, "y_max": 5, "y_axis_height": 10,
        "axes_color": BLACK,
        "graph_origin": ORIGIN,
    }
    def setup(self):
        GraphScene.setup(self)
        MovingCameraScene.setup(self)
    
    def construct(self):
        self.setup_axes()
        self.play(FadeOut(self.axes))
        self.wait()
        planegrid = VGroup()
        for i in range(-10, 10 + 1):
            planegrid.add(Line(5 * UP + i * RIGHT, 5 * DOWN + i * RIGHT))
        for i in range(-5, 5 + 1):
            planegrid.add(Line(i * UP + 10 * RIGHT, i * UP + 10 * LEFT))
        planegrid.set_color(GREY)
        planegrid.fade(0.875)
        self.add(planegrid)
        self.wait()
        powpt = TextMobject("$\\underline{\\text{Power of a point}}$").scale(1.25)
        powpt.add_background_rectangle()
        self.play(Write(powpt))
        self.play(powpt.shift, 3.25 * UP)
        cdot = Dot().fade(1)
        crad = ValueTracker(2)
        pdot = Dot(1.17 * RIGHT).rotate_about_origin(234 * DEGREES)
        self.add_foreground_mobject(pdot)
        plabel = TexMobject("P").scale(0.75).next_to(pdot, direction = DOWN, buff = SMALL_BUFF)
        plabel.add_background_rectangle()
        plabel.add_updater(lambda m: m.next_to(pdot, direction = DOWN, buff = SMALL_BUFF))
        circ = Circle(arc_center = cdot.get_center(), radius = crad.get_value())
        circlabel = TexMobject("C").scale(0.75).move_to(circ.get_center() + crad.get_value() * np.array([0.8, -0.6, 0]) + MED_LARGE_BUFF * np.array([0.8, 0.6, 0]))
        circlabel.add_background_rectangle()
        circlabel.add_updater(lambda m: m.move_to(circ.get_center() + crad.get_value() * np.array([0.8, -0.6, 0]) + MED_LARGE_BUFF * np.array([0.8, 0.6, 0])))
        self.play(Write(circ), Write(pdot), Write(plabel), Write(circlabel))
        self.wait()
        thet = ValueTracker(0.125)
        def intercreate(obj, t):
            interline = Line(pdot.get_center() + 10 * (obj.point_from_proportion(t) - pdot.get_center()), pdot.get_center() - 10 * (obj.point_from_proportion(thet.get_value()) - pdot.get_center())).set_color(BLUE)
            interpots = Intersection(interline, obj)
            alabel = TexMobject("A").scale(0.75).next_to(interpots[0], direction = DOWN, buff = SMALL_BUFF).add_background_rectangle()
            blabel = TexMobject("B").scale(0.75).next_to(interpots[1], direction = DOWN, buff = SMALL_BUFF).add_background_rectangle()
            return VGroup(interline, interpots, alabel, blabel)
        def interupd(obj):
            tempt = thet.get_value()
            tcirc = Circle(arc_center = cdot.get_center(), radius = crad.get_value())
            tempgp = intercreate(tcirc, tempt)
            obj.become(VGroup(tcirc, tempgp))
        intergp = intercreate(circ, thet.get_value())
        self.play(Write(intergp))
        self.wait()
        pwrpt = VGroup(
            TexMobject("\\mathcal{P}(P,C)"),
            TexMobject("=PA \\cdot PB"),
            TexMobject("=PT^2=PT'^2"),
            TexMobject("=PR \\cdot PS"),
        )
        for txt in pwrpt:
            txt.scale_in_place(0.75)
        pwrpt.arrange(DOWN)
        for txt in pwrpt[1:]:
            txt.align_to(pwrpt[0], LEFT)
            txt.add_background_rectangle()
        pwrpt[0].move_to(6 * LEFT + 1 * UP + pwrpt[0].get_center() - pwrpt[0].get_critical_point(LEFT))
        pwrpt[1].move_to(pwrpt[0].get_critical_point(RIGHT) + pwrpt[1].get_center() - pwrpt[1][0][0].get_critical_point(LEFT) + 0.15 * RIGHT)
        pwrpt[2].move_to(pwrpt[0].get_critical_point(RIGHT) + pwrpt[2].get_center() - pwrpt[2][0][0].get_critical_point(LEFT) + 0.15 * RIGHT)
        pwrpt[3].move_to(pwrpt[0].get_critical_point(RIGHT) + pwrpt[3].get_center() - pwrpt[3][0][0].get_critical_point(LEFT) + 2.1 * RIGHT)
        pwrpt[2].shift(0.625 * DOWN)
        self.play(Write(pwrpt[0]))
        self.wait()
        self.play(Write(pwrpt[1]))
        self.wait()
        intergp_circ = VGroup(circ, intergp)
        self.play(
            pdot.shift, 2 * LEFT,
            cdot.shift, 2 * RIGHT,
            crad.increment_value, -0.5,
            thet.set_value, 0.2,
            UpdateFromFunc(intergp_circ, interupd),
        )
        self.wait()
        sline = Line(pdot.get_center() + 10 * (circ.point_from_proportion(0.05) - pdot.get_center()), pdot.get_center() - 10 * (circ.point_from_proportion(0.05) - pdot.get_center())).set_color(BLUE)
        interspts = Intersection(sline, circ)
        rlabel = TexMobject("R").scale(0.75).next_to(interspts[0], direction = DOWN, buff = SMALL_BUFF).add_background_rectangle()
        slabel = TexMobject("S").scale(0.75).next_to(interspts[1], direction = DOWN, buff = SMALL_BUFF).add_background_rectangle()
        self.play(
            Write(sline),
            Write(interspts),
            Write(rlabel),
            Write(slabel)
        )
        self.wait()
        self.play(Write(pwrpt[3]))
        self.wait()
        tanline = CircleTangentLines(circ, pdot).set_color(BLUE)
        tlabel = TexMobject("T").scale(0.75).next_to(tanline[0].get_end(), direction = DOWN, buff = SMALL_BUFF).add_background_rectangle()
        tdlabel = TexMobject("T'").scale(0.75).next_to(tanline[1].get_end(), direction = UP, buff = SMALL_BUFF).add_background_rectangle()
        tangrp = VGroup(tanline, Dot(tanline[0].get_end()), Dot(tanline[1].get_end()), tlabel, tdlabel)
        self.play(Write(tangrp))
        self.wait()
        self.play(Write(pwrpt[2]))
        self.wait()
        plabel.clear_updaters()
        circlabel.clear_updaters()
        self.play(
            *[FadeOut(obj)for obj in [circ, tangrp, intergp, pwrpt, sline, interspts, rlabel, slabel, pdot, plabel, circlabel]]
        )
        self.wait(5)

class RdicalAxes(GraphScene, MovingCameraScene):
    CONFIG = {
        "x_min": -10, "x_max": 10, "x_axis_width": 20,
        "y_min": -5, "y_max": 5, "y_axis_height": 10,
        "axes_color": BLACK,
        "graph_origin": ORIGIN,
    }
    def setup(self):
        GraphScene.setup(self)
        MovingCameraScene.setup(self)
    
    def construct(self):
        self.setup_axes()
        self.play(FadeOut(self.axes))
        self.wait()
        planegrid = VGroup()
        for i in range(-10, 10 + 1):
            planegrid.add(Line(5 * UP + i * RIGHT, 5 * DOWN + i * RIGHT))
        for i in range(-5, 5 + 1):
            planegrid.add(Line(i * UP + 10 * RIGHT, i * UP + 10 * LEFT))
        planegrid.set_color(GREY)
        planegrid.fade(0.875)
        self.add(planegrid)
        powpt = TextMobject("$\\underline{\\text{Power of a point}}$").scale(1.25).shift(3.25 * UP)
        powpt_brect = BackgroundRectangle(powpt, color = BLACK, opacity = 0.75)
        powptgp = VGroup(powpt_brect, powpt)
        self.add_foreground_mobject(powptgp)
        self.add(powptgp)
        self.wait()
        cen1 = Dot(DOWN).fade(1)
        cen2 = Dot(DOWN).fade(1)
        cen1.shift(0.5 * LEFT)
        cen2.shift(2 * RIGHT)
        r1, r2 = ValueTracker(1.5), ValueTracker(2)
        c1 = Circle(arc_center = cen1.get_center(), radius = r1.get_value())
        c2 = Circle(arc_center = cen2.get_center(), radius = r2.get_value())
        self.play(Write(c1), Write(c2))
        self.wait()
        interpts = Intersection(c1, c2)
        alabel = TexMobject("A").scale(0.75).next_to(interpts[0], direction = RIGHT, buff = SMALL_BUFF).add_background_rectangle()
        blabel = TexMobject("B").scale(0.75).next_to(interpts[1], direction = RIGHT, buff = SMALL_BUFF).add_background_rectangle()
        ablabel = VGroup(alabel, blabel)
        radaxis = RadicalAxis(c1, c2)
        radaxis = Line(radaxis.get_start() + 10 * (radaxis.get_end() - radaxis.get_start()), radaxis.get_start() - 10 * (radaxis.get_end() - radaxis.get_start())).set_color(BLUE)
        self.play(Write(radaxis), Write(interpts), Write(ablabel))
        pdot = Intersection(radaxis, Line(ORIGIN, RIGHT)).shift(2.5 * UP)
        self.add_foreground_mobject(pdot)
        plabel = TexMobject("P").scale(0.75).next_to(pdot, direction = RIGHT, buff = SMALL_BUFF).add_background_rectangle()
        self.play(Write(pdot), Write(plabel))
        tanl1 = CircleTangentLines(c1, pdot)[0].set_color(BLUE)
        tanl2 = CircleTangentLines(c2, pdot)[1].set_color(BLUE)
        tlabel = TexMobject("T").scale(0.75).next_to(tanl1.get_end(), direction = LEFT, buff = MED_SMALL_BUFF).add_background_rectangle()
        tdlabel = TexMobject("T'").scale(0.75).next_to(tanl2.get_end(), direction = RIGHT, buff = MED_SMALL_BUFF).add_background_rectangle()
        tangrp = VGroup(tanl1, tanl2, Dot(tanl1.get_end()), Dot(tanl2.get_end()), tlabel, tdlabel)
        self.play(Write(tangrp))
        self.wait()
        pwrpt = TexMobject("PT^2=PA\\cdot PB = PT'^2").scale(0.75)
        pwrpt.move_to(6 * LEFT + 1 * UP + pwrpt.get_center() - pwrpt.get_critical_point(LEFT))
        self.play(Write(pwrpt))
        self.wait()
        eqtanl = TexMobject("\\implies PT = PT'").scale(0.75)
        eqtanl.move_to(pwrpt.get_center() + 1 * DOWN)
        self.play(Write(eqtanl))
        self.wait()
        self.play(
            *[FadeOut(obj) for obj in [eqtanl, pwrpt, tangrp, ablabel, pdot, plabel, interpts, radaxis]]
        )
        self.wait()
        c1.add_updater(lambda m: m.become(Circle(arc_center = cen1.get_center(), radius = r1.get_value())))
        self.play(cen1.shift, 2 * LEFT)
        rcircle = Circle(arc_center = 1.5 * DOWN + 0.5 * LEFT).set_color(GREEN)
        self.play(Write(rcircle))
        self.wait()
        intergp = VGroup(
            Intersection(rcircle, c1),
            Intersection(rcircle, c2),
        )
        self.add_foreground_mobject(intergp)
        labelgp = VGroup()
        alabel = TexMobject("A").scale(0.75).next_to(intergp[0][0], direction = LEFT, buff = SMALL_BUFF).add_background_rectangle()
        blabel = TexMobject("B").scale(0.75).next_to(intergp[0][1], direction = LEFT, buff = SMALL_BUFF).add_background_rectangle()
        rlabel = TexMobject("R").scale(0.75).next_to(intergp[1][0], direction = RIGHT, buff = SMALL_BUFF).add_background_rectangle()
        slabel = TexMobject("S").scale(0.75).next_to(intergp[1][1], direction = RIGHT, buff = SMALL_BUFF).add_background_rectangle()
        labelgp.add(alabel, blabel, rlabel, slabel)
        self.play(Write(intergp), Write(labelgp))
        self.wait()
        radgrp = VGroup()
        temp = RadicalAxis(rcircle, c1)
        temp = Line(temp.get_start() + 10 * (temp.get_end() - temp.get_start()), temp.get_start() - 10 * (temp.get_end() - temp.get_start())).set_color(BLUE)
        radgrp.add(temp)
        temp = RadicalAxis(rcircle, c2)
        temp = Line(temp.get_start() + 10 * (temp.get_end() - temp.get_start()), temp.get_start() - 10 * (temp.get_end() - temp.get_start())).set_color(BLUE)
        radgrp.add(temp)
        self.play(Write(radgrp))
        self.wait()
        radpt = Intersection(radgrp[0], radgrp[1])
        self.add_foreground_mobject(radpt)
        plabel = TexMobject("P").scale(0.75).next_to(radpt, direction = RIGHT, buff = SMALL_BUFF).add_background_rectangle()
        self.play(Write(radpt), Write(plabel))
        self.wait()
        tangrp = VGroup(
            CircleTangentLines(c1, radpt)[0],
            CircleTangentLines(c2, radpt)[1],
        ).set_color(BLUE)
        tangrp.add(Dot(tangrp[0].get_end()))
        tangrp.add(Dot(tangrp[1].get_end()))
        tlabel = TexMobject("T_1").scale(0.75).next_to(tangrp[-2].get_center(), direction = UP, buff = MED_SMALL_BUFF).add_background_rectangle()
        tdlabel = TexMobject("T_2").scale(0.75).next_to(tangrp[-1].get_center(), direction = UP, buff = MED_SMALL_BUFF).add_background_rectangle()
        self.play(Write(tangrp), Write(tlabel), Write(tdlabel))
        self.wait()
        pwrtxt = TexMobject("PT_1^2=PA \\cdot PB = PS \\cdot PR= PT_2^2").scale(0.75)
        pwrtxt.shift(4 * DOWN)
        pwrtxt[0][0:4].set_color(RED)
        pwrtxt[0][5:10].set_color(GREEN)
        pwrtxt[0][11:16].set_color(GREEN)
        pwrtxt[0][17:].set_color(RED)
        pwrtxt.add_background_rectangle()
        self.play(
            Write(pwrtxt),
            self.camera_frame.shift, DOWN,
        )
        self.wait()
        cengrp = VGroup(Dot(c1.get_center()), Dot(c2.get_center()))
        cenline = Line(cengrp[0].get_center(), cengrp[1].get_center())
        self.play(
            FadeOutAndShiftDown(pwrtxt),
            *[FadeOut(obj) for obj in [tlabel, tdlabel, tangrp, radgrp, labelgp, rcircle, intergp]],
            Write(cenline),
            Write(cengrp),
            self.camera_frame.shift, UP,
        )
        self.wait()
        #radline = Line(radpt.get_center(), cenline.get_start())
        radline = PerpendicularLine(radpt[0], cenline)
        radline = Line(radline.get_start() + 10 * (radline.get_end() - radline.get_start()), radline.get_start() - 10 * (radline.get_end() - radline.get_start())).set_color(BLUE)
        perp = Angle(cenline.get_start(), PerpendicularFoot(radpt[0], cenline).get_center(), radpt.get_center(), color = YELLOW, radius = 0.25)
        self.play(Write(radline), Write(perp))
        self.wait()
        c1.clear_updaters()
        self.play(FadeOut(perp), FadeOut(radpt), FadeOut(plabel))
        self.wait()
        randptpos = ValueTracker(0.45)
        def tangenerate(pos, radline0, obj):
            randpt = Dot(radline0.point_from_proportion(pos))
            #self.add_foreground_mobject(randpt)
            tanlen = Line(obj.get_center(), randpt.get_center()).get_length()
            if tanlen >= r1.get_value():
                if tanlen == 0:
                    tanlen += 0.001
                tanline = CircleTangentLines(obj, randpt)[0].set_color(BLUE)
            else:
                temp = Intersection(PerpendicularLine(randpt, Line(randpt.get_center(), obj.get_center())), obj)
                tanline = Line(randpt.get_center(), temp[0].get_center()).set_color(BLUE)
            tangrp = VGroup(tanline, Circle(arc_center = tanline.get_start(), radius = tanline.get_length()).set_color(PINK), Dot(tanline.get_end()))
            #self.add_foreground_mobject(tangrp[2])
            tangrp.add(Circle_Inversion(obj, about_point = tanline.get_start(), inversion_radius = tanline.get_length()))
            tangrp.add(Circle_Inversion(c2, about_point = tanline.get_start(), inversion_radius = tanline.get_length()))
            tanlen = Line(c2.get_center(), randpt.get_center()).get_length()
            if tanlen >= r2.get_value():
                if tanlen == 0:
                    tanlen += 0.001
                tanline = CircleTangentLines(c2, randpt)[1].set_color(BLUE)
            else:
                temp = Intersection(PerpendicularLine(randpt, Line(randpt.get_center(), c2.get_center())), c2)
                tanline = Line(randpt.get_center(), temp[0].get_center()).set_color(BLUE)
            tangrp.add(tanline, Dot(tanline.get_end()), randpt)
            return tangrp
        def tanupdate(obj):
            tpos = randptpos.get_value()
            tc1 = Circle(arc_center = cen1.get_center(), radius = r1.get_value())
            tcengrp = VGroup(Dot(tc1.get_center()), Dot(c2.get_center()))
            tcenline = Line(tcengrp[0].get_center(), tcengrp[1].get_center())
            tradline = RadicalAxis(tc1, c2)
            tradline = Line(tradline.get_start() + 10 * (tradline.get_end() - tradline.get_start()), tradline.get_start() - 10 * (tradline.get_end() - tradline.get_start())).set_color(BLUE)
            tempgrp = tangenerate(tpos, tradline, tc1)
            obj.become(VGroup(tc1, tcengrp, tcenline, tradline, tempgrp))
        ttangrp = tangenerate(randptpos.get_value(), radline, c1)
        self.play(Write(ttangrp))
        tttangrp = VGroup(c1, cengrp, cenline, radline, ttangrp)
        tttangrp.add_updater(tanupdate)
        self.add(tttangrp)
        self.play(
            randptpos.increment_value, 0.2,
            #UpdateFromFunc(tttangrp, tanupdate),
            run_time = 3
        )
        self.wait()
        self.play(
            cen1.shift, 2 * RIGHT,
            run_time = 3
        )
        self.wait()
        self.play(
            randptpos.increment_value, -0.075,
            run_time = 5,
        )
        self.play(
            randptpos.increment_value, -0.15,
            run_time = 5,
        )
        tttangrp.clear_updaters()
        self.wait()
        self.play(
            *[FadeOut(obj) for obj in [tttangrp, c1, c2]]
        )
        self.wait()
        fnlstmt = TextMobject("A line in which every point has the same power with respect \\\\ to two given circle is called the Radical Axis of two circles")
        fnlstmt.shift(DOWN)
        fnlstmt[0][12:37].set_color(BLUE)
        fnlstmt[0][75:86].set_color(YELLOW)
        fnlstmt.add_background_rectangle()
        self.play(
            FadeInFrom(fnlstmt),
            self.camera_frame.shift, DOWN + DOWN / 8,
        )
        self.wait(5)

class RdicalAxesThreeCircles(GraphScene, MovingCameraScene):
    CONFIG = {
        "x_min": -10, "x_max": 10, "x_axis_width": 20,
        "y_min": -5, "y_max": 5, "y_axis_height": 10,
        "axes_color": BLACK,
        "graph_origin": ORIGIN,
    }
    def setup(self):
        GraphScene.setup(self)
        MovingCameraScene.setup(self)
    
    def construct(self):
        self.setup_axes()
        self.play(FadeOut(self.axes))
        self.wait()
        planegrid = VGroup()
        for i in range(-10, 10 + 1):
            planegrid.add(Line(5 * UP + i * RIGHT, 5 * DOWN + i * RIGHT))
        for i in range(-5, 5 + 1):
            planegrid.add(Line(i * UP + 10 * RIGHT, i * UP + 10 * LEFT))
        planegrid.set_color(GREY)
        planegrid.fade(0.875)
        self.add(planegrid)
        powpt = TextMobject("$\\underline{\\text{Power of a point}}$").scale(1.25).shift(3.25 * UP)
        powpt_brect = BackgroundRectangle(powpt, color = BLACK, opacity = 0.75)
        powptgp = VGroup(powpt_brect, powpt)
        self.add_foreground_mobject(powptgp)
        self.add(powptgp)
        self.wait()
        cen1 = Dot(2 * RIGHT + 1 * UP)
        cen2 = Dot(3 * LEFT + 0.5 * UP)
        cen3 = Dot(0 * RIGHT + 2.5 * DOWN)
        c1 = Circle(arc_center = cen1.get_center(), radius = 2)
        c2 = Circle(arc_center = cen2.get_center(), radius = 1.5)
        c3 = Circle(arc_center = cen3.get_center(), radius = 1)
        self.play(Write(c1), Write(c2), Write(c3))
        self.wait()
        rad12 = RadicalAxis(c1, c2)
        rad23 = RadicalAxis(c2, c3)
        rad31 = RadicalAxis(c3, c1)
        def extendline(obj):
            return Line(obj.get_start() + 10 * (obj.get_end() - obj.get_start()), obj.get_start() - 10 * (obj.get_end() - obj.get_start()))
        rad12 = extendline(rad12).set_color(BLUE)
        rad23 = extendline(rad23).set_color(BLUE)
        rad31 = extendline(rad31).set_color(BLUE)
        radcen = RadicalCenter(c1, c2, c3)
        for obj in [rad12, rad23, rad31, radcen]:
            self.play(Write(obj))
        self.add_foreground_mobject(radcen)
        self.wait()
        tanlines = VGroup()
        for obj in [c1, c2, c3]:
            temptans = CircleTangentLines(obj, radcen).set_color(GREEN)
            tempgp = VGroup()
            tempgp.add(temptans)
            tempgp.add(Dot(temptans[0].get_end()))
            tempgp.add(Dot(temptans[1].get_end()))
            tanlines.add(tempgp)
        for obj in tanlines:
            self.play(Write(obj))
        radcircle = Circle(arc_center = radcen.get_center(), radius = tanlines[0][0][0].get_length()).set_color(PURPLE)
        self.play(Write(radcircle))
        self.wait()
        self.play(
            *[FadeOut(obj) for obj in [rad12, rad23, rad31, tanlines]]
        )
        self.wait()
        invcircgp = VGroup(
            Circle_Inversion(c1, about_point = radcen.get_center(), inversion_radius = tanlines[0][0][0].get_length()),
            Circle_Inversion(c2, about_point = radcen.get_center(), inversion_radius = tanlines[0][0][0].get_length()),
            Circle_Inversion(c3, about_point = radcen.get_center(), inversion_radius = tanlines[0][0][0].get_length()),
        )
        self.play(Write(invcircgp))
        self.wait()
        sols = ApolloniusCircles(c1, c2, c3)
        self.play(Write(sols[0]))
        self.wait()
        self.play(Write(sols[1]))
        self.wait()
        tanptsgrp = VGroup()
        for obj in [c1, c2, c3]:
            tempgp = VGroup(Intersection(obj, sols[0]), Intersection(obj, sols[1]))
            #tempgp = VGroup(HomotheticCenter(obj, sols[0])[0], HomotheticCenter(obj, sols[1])[1])
            tanptsgrp.add(tempgp)
            self.play(Write(tempgp))
        self.add_foreground_mobject(tanptsgrp)
        self.add_foreground_mobject(powptgp)
        invlinegp = VGroup()
        for j in range(3):
            invlinegp.add(Line(radcen.get_center(), tanptsgrp[j][0].get_center()))
            invlinegp.add(Line(radcen.get_center(), tanptsgrp[j][1].get_center()))
        invlinegp.set_color(BLUE)
        self.play(Write(invlinegp))
        self.wait()
        fnlstmt = VGroup(
            TextMobject("Radical center is a point which has the same power with \\\\ respect to three given circles"),
            TextMobject("Any circle mutually tangent to three given circles has a \\\\ conjugate circle that is also mutually tangent to the circles"),
            TextMobject("Both these circles are inverses of each other under an \\\\ inversion centered at the Radical center"),
            TextMobject("The point of tangency of each of these circles with the given \\\\ circles are also inverses of each other under the same inversion")
        )
        #self.add_foreground_mobject(fnlstmt)
        fnlstmt[0][0][:13].set_color(BLUE)
        fnlstmt[0][0][32:41].set_color(BLUE)
        #fnlstmt[0][0][59:].set_color(BLUE)
        fnlstmt[1][0][:9].set_color(GREEN)
        fnlstmt[1][0][47:62].set_color(GREEN)
        fnlstmt[2][0][19:38].set_color(ORANGE)
        fnlstmt[3][0][3:18].set_color(YELLOW)
        fnlstmt[3][0][64:83].set_color(YELLOW)
        fnlstmt.arrange(DOWN, buff = 0.875)
        fnlstmt.shift(DOWN)
        self.wait()
        poletitle = TextMobject("$\\underline{\\text{Poles and Polars}}$").scale(1.25).shift(3.25 * UP)
        poletitle.add_background_rectangle()
        self.play(
            *[FadeOut(obj) for obj in [tanptsgrp, invlinegp, invcircgp, radcen, sols[0], sols[1], c1, c2, c3, radcircle]],
            self.camera_frame.shift, DOWN + DOWN / 8,
        )
        self.wait()
        for txt in fnlstmt:
            txt.add_background_rectangle()
            self.play(FadeInFrom(txt))
        self.wait()
        self.play(
            FadeOutAndShift(fnlstmt),
            Transform(powpt, poletitle),
            self.camera_frame.shift, UP + UP / 8,
        )
        self.wait(5)

class PolesandPolars(GraphScene, MovingCameraScene):
    CONFIG = {
        "x_min": -10, "x_max": 10, "x_axis_width": 20,
        "y_min": -5, "y_max": 5, "y_axis_height": 10,
        "axes_color": BLACK,
        "graph_origin": ORIGIN,
    }
    def setup(self):
        GraphScene.setup(self)
        MovingCameraScene.setup(self)
    
    def construct(self):
        self.setup_axes()
        self.play(FadeOut(self.axes))
        self.wait()
        planegrid = VGroup()
        for i in range(-10, 10 + 1):
            planegrid.add(Line(5 * UP + i * RIGHT, 5 * DOWN + i * RIGHT))
        for i in range(-5, 5 + 1):
            planegrid.add(Line(i * UP + 10 * RIGHT, i * UP + 10 * LEFT))
        planegrid.set_color(GREY)
        planegrid.fade(0.875)
        self.add(planegrid)
        powpt = TextMobject("$\\underline{\\text{Poles and Polars}}$").scale(1.25).shift(3.25 * UP)
        powpt_brect = BackgroundRectangle(powpt, color = BLACK, opacity = 0.75)
        powptgp = VGroup(powpt_brect, powpt)
        self.add_foreground_mobject(powptgp)
        self.add(powptgp)
        self.wait()
        def extendline(obj):
            return Line(obj.get_start() + 10 * (obj.get_end() - obj.get_start()), obj.get_start() - 10 * (obj.get_end() - obj.get_start()))
        c1 = Circle(radius = 1.5)
        self.play(Write(c1))
        self.wait()
        def polecreator(obj):
            invpt = Circle_Inversion(obj, about_point = c1.get_center(), inversion_radius = Line(c1.get_center(), c1.get_start()).get_length())
            temppoleline = PerpendicularLine(invpt, Line(c1.get_center(), obj.get_center()))
            temppoleline = Line(temppoleline.get_start() + 10 * (temppoleline.get_end() - temppoleline.get_start()), temppoleline.get_start() - 10 * (temppoleline.get_end() - temppoleline.get_start())).set_color(BLUE)
            connecline = DashedLine(obj.get_center(), invpt.get_center()).set_color(BLUE)
            ang = Angle(temppoleline.get_end(), invpt.get_center(), obj.get_center(), radius = 0.25, color = YELLOW)
            return VGroup(ang, connecline, temppoleline)
        ppoint = Dot(4 * RIGHT + 2 * DOWN).set_color(PURPLE)
        plinegrp = polecreator(ppoint)
        for obj in plinegrp[1:]:
            obj.set_color(PURPLE)
        tanglines = CircleTangentLines(c1, ppoint).set_color(GREEN)
        tangp = VGroup()
        for obj in tanglines:
            tangp.add(extendline(obj).set_color(GREEN))
            tangp.add(Dot(obj.get_end()))
        self.play(Write(plinegrp[-1]))
        self.wait()
        self.play(Write(tangp))
        self.wait()
        self.play(Write(ppoint))
        self.wait()
        plabel = TexMobject("p").scale(0.75).move_to(3 * DOWN + 9 * LEFT / 8).add_background_rectangle()
        Plabel = TexMobject("P").scale(0.75).next_to(ppoint, direction = RIGHT, buff = SMALL_BUFF).add_background_rectangle()
        self.play(FadeOut(tangp), Write(plabel), Write(Plabel))
        self.wait()
        qpoint = Dot(5 * RIGHT + UP).set_color(BLUE)
        Qlabel = TexMobject("Q").scale(0.75).next_to(qpoint, direction = RIGHT, buff = SMALL_BUFF).add_background_rectangle()
        Qlabel.add_updater(lambda m: m.next_to(qpoint, direction = DOWN, buff = SMALL_BUFF))
        qlinegrp = polecreator(qpoint)
        qlabel = TexMobject("q").scale(0.75).next_to(qpoint, direction = RIGHT, buff = SMALL_BUFF).add_background_rectangle()
        qlabel.add_updater(lambda m: m.next_to(qlinegrp[-1].get_center() + LEFT / 4, direction = DOWN, buff = SMALL_BUFF))
        self.play(Write(qpoint), Write(Qlabel))
        self.wait()
        tanqlines = CircleTangentLines(c1, qpoint).set_color(GREEN)
        tanqlines.add(Dot(tanqlines[0].get_end()))
        tanqlines.add(Dot(tanqlines[1].get_end()))
        self.play(Write(tanqlines))
        self.wait()
        self.play(Write(qlinegrp[-1]), Write(qlabel))
        self.wait()
        #self.play(FadeOut(tanqlines))
        #self.wait()
        def qlineupd(obj):
            temp = polecreator(qpoint)[-1]
            obj.become(temp)
        def tanqupd(obj):
            ttanqlines = CircleTangentLines(c1, qpoint).set_color(GREEN)
            ttanqlines.add(Dot(ttanqlines[0].get_end()))
            ttanqlines.add(Dot(ttanqlines[1].get_end()))
            obj.become(ttanqlines)
        qlinegrp[-1].add_updater(qlineupd)
        self.add(qlinegrp[-1])
        tanqlines.add_updater(tanqupd)
        self.add(tanqlines)
        unitvec = (plinegrp[-1].get_end() - plinegrp[-1].get_start()) / plinegrp[-1].get_length()
        temppt = plinegrp[-1].get_start() + 7.75 * unitvec
        self.play(
            qpoint.move_to, temppt,
            run_time = 2
        )
        self.wait()
        self.play(
            qpoint.shift, 4.5 * unitvec,
            run_time = 5
        )
        self.wait()
        fpt = qpoint.get_center()
        path = VMobject()
        path.set_points_smoothly([fpt, [0, 2.5, 0], [-1, 0, 0], [0, -1, 0], [3, 0, 0], fpt])
        #self.add(path)
        self.play(
            MoveAlongPath(qpoint, path),
            run_time = 15,
            rate_func = linear
        )
        qlinegrp[-1].clear_updaters()
        tanqlines.clear_updaters()
        self.play(
            VGroup(tanqlines, c1, plabel, qlabel, Plabel, Qlabel, ppoint, qpoint, qlinegrp[-1], plinegrp[-1]).fade, 0.875,
        )
        fnlstmt = VGroup(
            TextMobject("If the pole $P$ of a line $p$ lies on another line $q$, \\\\ then the pole $Q$ of the line $q$ lies on line $p$"),
            TextMobject("If the polar $p$ of a point $P$ passes through another point $Q$, \\\\ then the polar $q$ of point $Q$ passes through point $P$")
        ).arrange(DOWN, buff = 1)
        fnlstmt[0][0][5:10].set_color(RED)
        fnlstmt[0][0][13:18].set_color(RED)
        fnlstmt[0][0][31:36].set_color(BLUE)
        fnlstmt[0][0][44:49].set_color(BLUE)
        fnlstmt[0][0][54:59].set_color(BLUE)
        fnlstmt[0][0][65:70].set_color(RED)
        fnlstmt[1][0][5:11].set_color(RED)
        fnlstmt[1][0][14:20].set_color(RED)
        fnlstmt[1][0][40:46].set_color(BLUE)
        fnlstmt[1][0][54:60].set_color(BLUE)
        fnlstmt[1][0][62:68].set_color(BLUE)
        fnlstmt[1][0][81:87].set_color(RED)
        for txt in fnlstmt:
            txt.add_background_rectangle()
            self.play(Write(txt))
        self.wait(5)

class HomoCenters(GraphScene, MovingCameraScene):
    CONFIG = {
        "x_min": -10, "x_max": 10, "x_axis_width": 20,
        "y_min": -5, "y_max": 5, "y_axis_height": 10,
        "axes_color": BLACK,
        "graph_origin": ORIGIN,
    }
    def setup(self):
        GraphScene.setup(self)
        MovingCameraScene.setup(self)
    
    def construct(self):
        self.setup_axes()
        self.play(FadeOut(self.axes))
        self.wait()
        planegrid = VGroup()
        for i in range(-10, 10 + 1):
            planegrid.add(Line(5 * UP + i * RIGHT, 5 * DOWN + i * RIGHT))
        for i in range(-5, 5 + 1):
            planegrid.add(Line(i * UP + 10 * RIGHT, i * UP + 10 * LEFT))
        planegrid.set_color(GREY)
        planegrid.fade(0.875)
        self.add(planegrid)
        c1 = Circle(arc_center = 1 * LEFT)
        c2 = Circle(arc_center = 4 * RIGHT, radius = 2)
        cengrp = VGroup(Dot(c1.get_center()), Dot(c2.get_center()))
        c1label = TexMobject("C_1").scale(0.75).next_to(cengrp[0], direction = DOWN, buff = SMALL_BUFF).add_background_rectangle()
        c2label = TexMobject("C_2").scale(0.75).next_to(cengrp[1], direction = DOWN, buff = SMALL_BUFF).add_background_rectangle()
        cenlabel = VGroup(c1label, c2label)
        self.play(Write(c1), Write(c2), Write(cengrp), Write(cenlabel))
        self.wait()
        h = HomotheticCenter(c1, c2)
        elabel = TexMobject("E").scale(0.75).next_to(h[1], direction = DOWN, buff = SMALL_BUFF).add_background_rectangle()
        ilabel = TexMobject("I").scale(0.75).next_to(h[0], direction = DOWN, buff = SMALL_BUFF).add_background_rectangle()
        ename = TextMobject("''External''", " Homothetic Center").scale(0.75).move_to(4 * LEFT + 3 * UP)
        iname = TextMobject("''Internal''", " Homothetic Center").scale(0.75).move_to(4 * LEFT + 3 * UP)
        ename.set_color_by_tex("External", YELLOW)
        iname.set_color_by_tex("Internal", YELLOW)
        ename.add_background_rectangle()
        iname.add_background_rectangle()
        exttans = VGroup()
        tanpts = VGroup()
        tgrp = VGroup()
        for obj in CircleTangentLines(c2, h[1]):
            temp = Line(obj.get_start() + 20 * (obj.get_end() - obj.get_start()) / obj.get_length(), obj.get_start() - 20 * (obj.get_end() - obj.get_start()) / obj.get_length()).set_color(GREEN)
            tgrp.add(Dot(obj.get_end()))
            exttans.add(temp)
        self.play(Write(exttans))
        self.wait()
        self.play(Write(h[1]), Write(elabel), Write(ename))
        self.add_foreground_mobject(h[1])
        tanpts.add(tgrp)
        tgrp = VGroup()
        inttans = VGroup()
        for obj in CircleTangentLines(c2, h[0]):
            temp = Line(obj.get_start() + 20 * (obj.get_end() - obj.get_start()) / obj.get_length(), obj.get_start() - 20 * (obj.get_end() - obj.get_start()) / obj.get_length()).set_color(GREEN)
            tgrp.add(Dot(obj.get_end()))
            inttans.add(temp)
        tanpts.add(tgrp)
        tgrp = VGroup()
        for obj in CircleTangentLines(c1, h[1]):
            tgrp.add(Dot(obj.get_end()))
        tanpts.add(tgrp)
        tgrp = VGroup()
        for obj in CircleTangentLines(c1, h[0]):
            tgrp.add(Dot(obj.get_end()))
        tanpts.add(tgrp)
        self.play(Write(tanpts[0]), Write(tanpts[2]))
        self.wait()
        thet = ValueTracker(1 / 8)
        def tancriccreator(t, p):
            if t %  (1 / 2) == 0:
                t -= 1 / 2048
            cpt = Line(c2.get_center(), c2.get_start())
            cpt.rotate(t * 2 * PI, about_point = c2.get_center())
            cpt = Line(h[p].get_center(), cpt.get_end())
            cpt = Intersection(cpt, c2).set_color(GREEN)
            #self.play(Write(cpt))
            ocpt = Line(c2.get_center(), c2.get_start())
            ocpt.rotate(t * 2 * PI, about_point = c2.get_center())
            ocpt = Line(h[p].get_center(), ocpt.get_end())
            ocpt = Intersection(ocpt, c1).set_color(GREEN)
            #self.play(Write(cpt))
            fcen = Intersection(Line(c2.get_center(), cpt[0].get_center()), Line(c1.get_center(), ocpt[(0 + p) % 2].get_center())).set_color(PINK)
            scen = Intersection(Line(c2.get_center(), cpt[1].get_center()), Line(c1.get_center(), ocpt[(1 + p) % 2].get_center())).set_color(PINK)
            #self.play(Write(fcen), Write(scen), Write(cpt), Write(ocpt))
            resgrp = VGroup()
            resgrp.add(Circle(arc_center = fcen.get_center(), radius = Line(fcen.get_center(), cpt[0].get_center()).get_length()).set_color(GREEN))
            resgrp.add(Circle(arc_center = scen.get_center(), radius = Line(scen.get_center(), cpt[1].get_center()).get_length()).set_color(GREEN))
            return resgrp
        def circlefamily(p, dt = 0.125):
            resgrp = VGroup()
            for t in np.arange(0, 1, dt):
                resgrp.add(tancriccreator(t, p))
            return resgrp
        def extupd(obj):
            temp = tancriccreator(thet.get_value(), 1)
            obj.become(temp)
        def intupd(obj):
            temp = tancriccreator(thet.get_value(), 0)
            obj.become(temp)
        c1cpy = c1.copy()
        self.add(c1cpy)
        self.play(
            c1cpy.move_to, c2.get_center(),
            c1cpy.scale, 2,
            FadeOut(exttans),
            FadeOut(tanpts[0]),
            FadeOut(tanpts[2]),
            #c1cpy.fade, 0.5,
            run_time = 2
        )
        self.remove(c1cpy)
        self.wait()
        def tancriccreator2(t, p):
            if t %  (1 / 2) == 0:
                t -= 1 / 2048
            cpt = Line(c2.get_center(), c2.get_start())
            cpt.rotate(t * 2 * PI, about_point = c2.get_center())
            cpt = Line(h[p].get_center(), cpt.get_end())
            cptcpy = cpt.copy()
            cpt = Intersection(cpt, c2)
            #self.play(Write(cpt))
            ocpt = Line(c2.get_center(), c2.get_start())
            ocpt.rotate(t * 2 * PI, about_point = c2.get_center())
            ocpt = Line(h[p].get_center(), ocpt.get_end())
            ocpt = Intersection(ocpt, c1)
            #self.play(Write(cpt))
            fcen = Intersection(Line(c2.get_center(), cpt[0].get_center()), Line(c1.get_center(), ocpt[(0 + p) % 2].get_center())).set_color(PINK)
            scen = Intersection(Line(c2.get_center(), cpt[1].get_center()), Line(c1.get_center(), ocpt[(1 + p) % 2].get_center())).set_color(PINK)
            #self.play(Write(fcen), Write(scen), Write(cpt), Write(ocpt))
            resgrp = VGroup()
            resgrp.add(Circle(arc_center = fcen.get_center(), radius = Line(fcen.get_center(), cpt[0].get_center()).get_length()).set_color(GREEN))
            resgrp.add(Circle(arc_center = scen.get_center(), radius = Line(scen.get_center(), cpt[1].get_center()).get_length()).set_color(GREEN))
            cptcpy = Line(cptcpy.get_start(), cptcpy.get_start() + (cptcpy.get_end() - cptcpy.get_start()) / cptcpy.get_length())
            cptcpy = Line(cptcpy.get_start() - 10 * (cptcpy.get_end() - cptcpy.get_start()), cptcpy.get_start() + 20 * (cptcpy.get_end() - cptcpy.get_start())).set_color(BLUE)
            resgrp.add(cptcpy, cpt, ocpt)
            return resgrp
        def extupd2(obj):
            temp = tancriccreator2(thet.get_value(), 1)
            obj.become(temp)
        def intupd2(obj):
            temp = tancriccreator2(thet.get_value(), 0)
            obj.become(temp)
        examplines = tancriccreator2(thet.get_value(), 1)
        alabel = TexMobject("A").scale(0.75).next_to(examplines[-1][1], direction = LEFT, buff = SMALL_BUFF).add_background_rectangle()
        blabel = TexMobject("B").scale(0.75).next_to(examplines[-1][0].get_center(), direction = DR, buff = SMALL_BUFF).add_background_rectangle()
        adlabel = TexMobject("A'").scale(0.75).next_to(examplines[-2][1].get_center(), direction = DL, buff = SMALL_BUFF).add_background_rectangle()
        bdlabel = TexMobject("B'").scale(0.75).next_to(examplines[-2][0], direction = RIGHT, buff = SMALL_BUFF).add_background_rectangle()
        labelgrp = VGroup(alabel, blabel, adlabel, bdlabel)
        self.play(Write(examplines[2:]), Write(labelgrp))
        #self.add_foreground_mobject(labelgrp)
        self.wait()
        trigrp = VGroup(
            Line(c1.get_center(), examplines[-1][1].get_center()),
            Line(c1.get_center(), examplines[-1][0].get_center()),
            Line(c2.get_center(), examplines[-2][1].get_center()),
            Line(c2.get_center(), examplines[-2][0].get_center()),
            Angle(c1.get_center(), examplines[-1][1].get_center(), examplines[-1][0].get_center(), radius = 0.375),
            Angle(c1.get_center(), examplines[-1][0].get_center(), examplines[-1][1].get_center(), radius = 0.375),
            Angle(c2.get_center(), examplines[-2][1].get_center(), examplines[-2][0].get_center(), radius = 0.375),
            Angle(c2.get_center(), examplines[-2][0].get_center(), examplines[-2][1].get_center(), radius = 0.375),
            TexMobject("C_1A=C_1B \\implies \\angle C_1AB=\\angle C_1BA").scale(0.625).move_to(c1.get_center() + 1.5 * DOWN).add_background_rectangle(),
            TexMobject("C_2A'=C_2B' \\implies \\angle C_2A'B'=\\angle C_2B'A'").scale(0.625).move_to(c2.get_center() + 2.5 * DOWN).add_background_rectangle(),
        )
        trigrp[:4].set_color(BLUE)
        self.play(Write(trigrp))
        self.wait()
        pintersec = Intersection(Line(c1.get_center(), examplines[-1][0].get_center()), Line(c2.get_center(), examplines[-2][1].get_center()))
        plabel = TexMobject("P").scale(0.75).next_to(pintersec, direction = UP, buff = SMALL_BUFF).add_background_rectangle()
        ptrigrp = VGroup(Line(examplines[-1][0].get_center(), pintersec.get_center()), Line(examplines[-2][1].get_center(), pintersec.get_center()))
        ptrigrp.set_color(BLUE)
        ptrigrp.add(pintersec, plabel)
        self.play(Write(ptrigrp))
        self.wait()
        ptrigrp.add(
            Angle(examplines[-2][1].get_center(), examplines[-1][0].get_center(), pintersec.get_center(), radius = 0.375),
            Angle(examplines[-1][0].get_center(), examplines[-2][1].get_center(), pintersec.get_center(), radius = 0.375),
        )
        self.play(Write(ptrigrp[-2:]))
        self.wait()
        self.play(Write(examplines[1]), FadeOut(trigrp[-1]), FadeOut(trigrp[-2]))
        self.wait()
        qintersec = Intersection(Line(examplines[-1][1].get_center(), c1.get_center()), Line(examplines[-2][0].get_center(), c2.get_center()))
        qlabel = TexMobject("Q").scale(0.75).next_to(qintersec, direction = DOWN, buff = SMALL_BUFF).add_background_rectangle()
        qtrigrp = VGroup(
            Line(qintersec.get_center(), c1.get_center()),
            Line(qintersec.get_center(), c2.get_center()),
        ).set_color(BLUE)
        qtrigrp.add(qintersec, qlabel)
        self.play(Write(qtrigrp))
        self.wait()
        self.play(Write(examplines[0]))
        self.wait()
        self.play(
            FadeOut(trigrp[:-2]),
            FadeOut(qtrigrp),
            FadeOut(ptrigrp),
            FadeOut(labelgrp)
        )
        self.wait()
        self.play(
            thet.increment_value, 4 / 8 - 1 / 8,
            UpdateFromFunc(examplines, extupd2),
            run_time = 9,
            rate_func = there_and_back
        )
        self.wait()
        self.play(FadeOut(elabel), FadeOut(h[1]), FadeOut(examplines))
        self.play(Write(inttans), Write(tanpts[1]), Write(tanpts[3]))
        self.play(Write(h[0]), Write(ilabel), Transform(ename, iname))
        self.wait()
        c1cpy = c1.copy()
        self.add(c1cpy)
        scalefac = ValueTracker(0)
        c1cpy.add_updater(lambda m: m.become(Circle(arc_center = c1.get_center() + scalefac.get_value() * (c2.get_center() - c1.get_center()), radius = abs(3 * scalefac.get_value() - 1))))
        self.play(
            scalefac.set_value, 1,
            run_time = 3
        )
        c1cpy.clear_updaters()
        self.remove(c1cpy)
        self.wait()
        examplines = tancriccreator2(thet.get_value(), 0)
        self.play(
            Write(examplines),
            *[FadeOut(obj) for obj in [inttans, tanpts[1], tanpts[3]]]
        )
        self.wait()
        self.play(
            thet.increment_value, 4 / 8 - 1 / 8,
            UpdateFromFunc(examplines, intupd2),
            run_time = 9,
            rate_func = there_and_back
        )
        self.wait()
        self.play(FadeOut(examplines), FadeOut(h[0]), FadeOut(ilabel), FadeOut(ename))
        self.wait()
        cfam = circlefamily(1, dt = 1 / 62)
        for obj in cfam:
            self.add(obj)
            self.wait(1 / 8)
        self.wait()
        self.play(FadeOut(cfam))
        self.wait()
        cfam = circlefamily(0, dt = 1 / 64)
        for obj in cfam:
            self.add(obj)
            self.wait(1 / 8)
        self.wait()
        self.play(FadeOut(cfam))
        self.wait()
        '''self.play(
            thet.set_value, 0.5,
            UpdateFromFunc(sols, extupd),
            run_time = 6
        )
        self.wait()
        thet.set_value(0)
        self.play(Transform(sols, tancriccreator(thet.get_value(), 0)))
        self.wait()
        self.play(
            thet.set_value, 0.5,
            UpdateFromFunc(sols, intupd),
            run_time = 6
        )'''
        self.wait(5)

class SamePowerInFamily(GraphScene, MovingCameraScene):
    CONFIG = {
        "x_min": -10, "x_max": 10, "x_axis_width": 20,
        "y_min": -5, "y_max": 5, "y_axis_height": 10,
        "axes_color": BLACK,
        "graph_origin": ORIGIN,
    }
    def setup(self):
        GraphScene.setup(self)
        MovingCameraScene.setup(self)
    
    def construct(self):
        self.setup_axes()
        self.play(FadeOut(self.axes))
        self.wait()
        planegrid = VGroup()
        for i in range(-10, 10 + 1):
            planegrid.add(Line(5 * UP + i * RIGHT, 5 * DOWN + i * RIGHT))
        for i in range(-5, 5 + 1):
            planegrid.add(Line(i * UP + 10 * RIGHT, i * UP + 10 * LEFT))
        planegrid.set_color(GREY)
        planegrid.fade(0.875)
        self.add(planegrid)
        c1 = Circle(arc_center = 2.5 * LEFT)
        c2 = Circle(arc_center = 2.5 * RIGHT, radius = 2)
        cengrp = VGroup(Dot(c1.get_center()), Dot(c2.get_center()))
        c1label = TexMobject("C_1").scale(0.75).next_to(cengrp[0], direction = DOWN, buff = SMALL_BUFF).add_background_rectangle()
        c2label = TexMobject("C_2").scale(0.75).next_to(cengrp[1], direction = DOWN, buff = SMALL_BUFF).add_background_rectangle()
        cenlabel = VGroup(c1label, c2label)
        self.add(c1, c2, cengrp, cenlabel)
        h = HomotheticCenter(c1, c2)
        elabel = TexMobject("E").scale(0.75).next_to(h[1], direction = DOWN, buff = SMALL_BUFF).add_background_rectangle()
        ilabel = TexMobject("I").scale(0.75).next_to(h[0], direction = DOWN, buff = SMALL_BUFF).add_background_rectangle()
        exttans = VGroup()
        tanpts = VGroup()
        tgrp = VGroup()
        for obj in CircleTangentLines(c2, h[1]):
            temp = Line(obj.get_start() + 20 * (obj.get_end() - obj.get_start()) / obj.get_length(), obj.get_start() - 20 * (obj.get_end() - obj.get_start()) / obj.get_length()).set_color(GREEN)
            tgrp.add(Dot(obj.get_end()))
            exttans.add(temp)
        tanpts.add(tgrp)
        tgrp = VGroup()
        inttans = VGroup()
        for obj in CircleTangentLines(c2, h[0]):
            temp = Line(obj.get_start() + 20 * (obj.get_end() - obj.get_start()) / obj.get_length(), obj.get_start() - 20 * (obj.get_end() - obj.get_start()) / obj.get_length()).set_color(GREEN)
            tgrp.add(Dot(obj.get_end()))
            inttans.add(temp)
        tanpts.add(tgrp)
        tgrp = VGroup()
        for obj in CircleTangentLines(c1, h[1]):
            tgrp.add(Dot(obj.get_end()))
        tanpts.add(tgrp)
        tgrp = VGroup()
        for obj in CircleTangentLines(c1, h[0]):
            tgrp.add(Dot(obj.get_end()))
        tanpts.add(tgrp)
        self.add(h[0], ilabel)
        self.add_foreground_mobject(h[0])
        self.wait()
        def tancriccreator(t, p):
            if t %  (1 / 2) == 0:
                t -= 1 / 2048
            cpt = Line(c2.get_center(), c2.get_start())
            cpt.rotate(t * 2 * PI, about_point = c2.get_center())
            cpt = Line(h[p].get_center(), cpt.get_end())
            cpt = Intersection(cpt, c2).set_color(GREEN)
            #self.play(Write(cpt))
            ocpt = Line(c2.get_center(), c2.get_start())
            ocpt.rotate(t * 2 * PI, about_point = c2.get_center())
            ocpt = Line(h[p].get_center(), ocpt.get_end())
            ocpt = Intersection(ocpt, c1).set_color(GREEN)
            #self.play(Write(cpt))
            fcen = Intersection(Line(c2.get_center(), cpt[0].get_center()), Line(c1.get_center(), ocpt[(0 + p) % 2].get_center())).set_color(PINK)
            scen = Intersection(Line(c2.get_center(), cpt[1].get_center()), Line(c1.get_center(), ocpt[(1 + p) % 2].get_center())).set_color(PINK)
            #self.play(Write(fcen), Write(scen), Write(cpt), Write(ocpt))
            resgrp = VGroup()
            resgrp.add(Circle(arc_center = fcen.get_center(), radius = Line(fcen.get_center(), cpt[0].get_center()).get_length()).set_color(GREEN))
            resgrp.add(Circle(arc_center = scen.get_center(), radius = Line(scen.get_center(), cpt[1].get_center()).get_length()).set_color(GREEN))
            return resgrp
        def circlefamily(p, dt = 0.125):
            resgrp = VGroup()
            for t in np.arange(0, 1, dt):
                resgrp.add(tancriccreator(t, p))
            return resgrp
        def extupd(obj):
            temp = tancriccreator(thet.get_value(), 1)
            obj.become(temp)
        def intupd(obj):
            temp = tancriccreator(thet.get_value(), 0)
            obj.become(temp)
        def tancriccreator2(t, p):
            if t %  (1 / 2) == 0:
                t -= 1 / 2048
            cpt = Line(c2.get_center(), c2.get_start())
            cpt.rotate(t * 2 * PI, about_point = c2.get_center())
            cpt = Line(h[p].get_center(), cpt.get_end())
            cptcpy = cpt.copy()
            cpt = Intersection(cpt, c2)
            #self.play(Write(cpt))
            ocpt = Line(c2.get_center(), c2.get_start())
            ocpt.rotate(t * 2 * PI, about_point = c2.get_center())
            ocpt = Line(h[p].get_center(), ocpt.get_end())
            ocpt = Intersection(ocpt, c1)
            #self.play(Write(cpt))
            fcen = Intersection(Line(c2.get_center(), cpt[0].get_center()), Line(c1.get_center(), ocpt[(0 + p) % 2].get_center())).set_color(PINK)
            scen = Intersection(Line(c2.get_center(), cpt[1].get_center()), Line(c1.get_center(), ocpt[(1 + p) % 2].get_center())).set_color(PINK)
            #self.play(Write(fcen), Write(scen), Write(cpt), Write(ocpt))
            resgrp = VGroup()
            resgrp.add(Circle(arc_center = fcen.get_center(), radius = Line(fcen.get_center(), cpt[0].get_center()).get_length()).set_color(GREEN))
            resgrp.add(Circle(arc_center = scen.get_center(), radius = Line(scen.get_center(), cpt[1].get_center()).get_length()).set_color(GREEN))
            cptcpy = Line(cptcpy.get_start(), cptcpy.get_start() + (cptcpy.get_end() - cptcpy.get_start()) / cptcpy.get_length())
            cptcpy = Line(cptcpy.get_start() - 10 * (cptcpy.get_end() - cptcpy.get_start()), cptcpy.get_start() + 20 * (cptcpy.get_end() - cptcpy.get_start())).set_color(BLUE)
            resgrp.add(cptcpy, cpt, ocpt)
            return resgrp
        def extupd2(obj):
            temp = tancriccreator2(thet.get_value(), 1)
            obj.become(temp)
        def intupd2(obj):
            temp = tancriccreator2(thet.get_value(), 0)
            obj.become(temp)
        randcirc = tancriccreator2(1 / 16, 0)
        plabel = TexMobject("P").scale(0.75).next_to(tanpts[3][0], direction = UP, buff = SMALL_BUFF).add_background_rectangle()
        qlabel = TexMobject("Q").scale(0.75).next_to(tanpts[1][0], direction = DOWN, buff = SMALL_BUFF).add_background_rectangle()
        alabel = TexMobject("A").scale(0.75).next_to(randcirc[-1][0], direction = DOWN, buff = SMALL_BUFF).add_background_rectangle()
        blabel = TexMobject("B").scale(0.75).next_to(randcirc[-2][0], direction = DOWN, buff = SMALL_BUFF).add_background_rectangle()
        adlabel = TexMobject("A'").scale(0.75).next_to(randcirc[-2][1], direction = DOWN, buff = SMALL_BUFF).add_background_rectangle()
        bdlabel = TexMobject("B'").scale(0.75).next_to(randcirc[-1][1], direction = DOWN, buff = SMALL_BUFF).add_background_rectangle()
        labelgrp = VGroup(alabel, blabel, adlabel, bdlabel, plabel, qlabel)
        self.play(Write(randcirc[0]), Write(randcirc[2:]), Write(labelgrp[:-2]))
        self.wait()
        self.play(Write(inttans[0]), Write(tanpts[1][0]), Write(tanpts[3][0]), Write(labelgrp[-2:]))
        self.wait()
        ipagrp = VGroup(
            Line(tanpts[3][0].get_center(), randcirc[-1][0].get_center()).set_color(BLUE),
            Angle(randcirc[-1][0].get_center(), tanpts[3][0].get_center(), h[0].get_center(), radius = 0.375, color = PURPLE),
            Line(tanpts[3][0].get_center(), randcirc[-1][1].get_center()).set_color(BLUE),
            Angle(tanpts[3][0].get_center(), randcirc[-1][1].get_center(), h[0].get_center(), radius = 0.375, color = PURPLE),
        )
        self.play(Write(ipagrp))
        self.wait()
        iqbgrp = VGroup(
            Line(randcirc[-2][1].get_center(), tanpts[1][0].get_center()).set_color(BLUE),
            Line(tanpts[1][0].get_center(), randcirc[-2][0].get_center()).set_color(BLUE),
            Angle(tanpts[1][0].get_center(), randcirc[-2][0].get_center(), randcirc[-2][1].get_center(), color = PURPLE, radius = 0.375)
        )
        self.play(Write(iqbgrp))
        self.wait()
        aqline = DashedLine(randcirc[-1][0].get_center(), tanpts[1][0].get_center()).set_color(YELLOW)
        self.play(Write(aqline))
        self.wait()
        '''indicangles = VGroup(
            VGroup(Line(randcirc[-1][0].get_center(), tanpts[3][0].get_center()), Line(tanpts[3][0].get_center(), tanpts[1][0].get_center()), Line(tanpts[1][0].get_center(), randcirc[-1][0].get_center())),
            VGroup(Line(randcirc[-1][0].get_center(), randcirc[-2][0].get_center()), Line(randcirc[-2][0].get_center(), tanpts[1][0].get_center()), Line(tanpts[1][0].get_center(), randcirc[-1][0].get_center()))
        ).set_color(YELLOW)'''
        indicangles = VGroup(
            VGroup(Line(randcirc[-1][0].get_center(), tanpts[3][0].get_center()), Line(tanpts[3][0].get_center(), tanpts[1][0].get_center())).set_color(GOLD),
            VGroup(Line(randcirc[-1][0].get_center(), randcirc[-2][0].get_center()), Line(randcirc[-2][0].get_center(), tanpts[1][0].get_center())).set_color(ORANGE)
        )
        self.play(Indicate(indicangles[0], scale_factor = 1.5))
        self.play(Indicate(indicangles[1], scale_factor = 1.5))
        #self.play(Indicate(indicangles, scale_factor = 1.5))
        self.play(FadeOut(indicangles), run_time = 0.01)
        #self.remove(indicangles)
        self.wait()
        fcir = ThreePointCircle(tanpts[3][0], randcirc[-1][0], tanpts[1][0]).set_color(GREY)
        #clabel = TexMobject("C").scale(0.75).next_to(fcir, direction = LEFT, buff = SMALL_BUFF).add_background_rectangle()
        clabel = TexMobject("C").scale(0.75).move_to(3 * UP + 2 * LEFT).add_background_rectangle()
        self.play(
            Write(fcir),
            Write(clabel),
            *[FadeOut(obj) for obj in [iqbgrp, ipagrp, aqline]]
        )
        self.wait()
        pwreqn = TexMobject("\\mathcal{P}(I,C) = IA \\cdot IB = IP \\cdot IQ").scale(0.75).add_background_rectangle()
        pwreqn.move_to(2.5 * DOWN + 6 * LEFT + pwreqn.get_center() - pwreqn.get_critical_point(LEFT))
        self.play(Write(pwreqn))
        self.wait()
        alabel.add_updater(lambda m: m.next_to(randcirc[-1][0], direction = DOWN, buff = SMALL_BUFF).add_background_rectangle())
        blabel.add_updater(lambda m: m.next_to(randcirc[-2][0], direction = DOWN, buff = SMALL_BUFF).add_background_rectangle())
        adlabel.add_updater(lambda m: m.next_to(randcirc[-2][1], direction = DOWN, buff = SMALL_BUFF).add_background_rectangle())
        bdlabel.add_updater(lambda m: m.next_to(randcirc[-1][1], direction = DOWN, buff = SMALL_BUFF).add_background_rectangle())
        thet = ValueTracker(1 / 16)
        rrandcirc = VGroup(randcirc[0], *randcirc[2:], fcir)
        def upd3(obj):
            trandcir = tancriccreator2(thet.get_value(), 0)
            tfcir = ThreePointCircle(tanpts[3][0], trandcir[-1][0], tanpts[1][0]).set_color(GREY)
            tgrp = VGroup(trandcir[0], *trandcir[2:], tfcir)
            obj.become(tgrp)
        self.play(FadeOut(clabel))
        self.play(
            thet.increment_value, 0.375,
            UpdateFromFunc(rrandcirc, upd3),
            run_time = 15,
            rate_func = there_and_back
        )
        alabel.clear_updaters()
        blabel.clear_updaters()
        adlabel.clear_updaters()
        bdlabel.clear_updaters()
        self.play(
            *[FadeOut(obj) for obj in [h[0], ilabel, pwreqn, tanpts[1][0], tanpts[3][0], rrandcirc, labelgrp, inttans[0]]]
        )
        flstmt = VGroup(
            TextMobject("Both the Homothetic centers have the same power with\\\\respect to all the circles in its own family."),
            TextMobject("In other words,"),
            TextMobject("Each Homothetic center lies on the radical axis of any two\\\\circles of its own family.")
        )
        flstmt.arrange(DOWN, buff = 1)
        flstmt[1].set_color(YELLOW)
        flstmt[0][0][7:24].set_color(BLUE)
        flstmt[0][0][31:40].set_color(BLUE)
        flstmt[0][0][68:80].set_color(BLUE)
        flstmt[2][0][4:20].set_color(BLUE)
        flstmt[2][0][29:40].set_color(BLUE)
        flstmt[2][0][57:69].set_color(BLUE)
        self.play(
            Write(flstmt),
            cengrp.fade, 1 - 1 / 8,
            cenlabel.fade, 1 - 1 / 8,
            c1.fade, 1 - 1 / 8,
            c2.fade, 1 - 1 / 8,
        )
        self.wait(5)

class HomotheticAxes(GraphScene, MovingCameraScene):
    CONFIG = {
        "x_min": -10, "x_max": 10, "x_axis_width": 20,
        "y_min": -5, "y_max": 5, "y_axis_height": 10,
        "axes_color": BLACK,
        "graph_origin": ORIGIN,
    }
    def setup(self):
        GraphScene.setup(self)
        MovingCameraScene.setup(self)
    
    def construct(self):
        self.setup_axes()
        self.play(FadeOut(self.axes))
        self.wait()
        planegrid = VGroup()
        for i in range(-10, 10 + 1):
            planegrid.add(Line(5 * UP + i * RIGHT, 5 * DOWN + i * RIGHT))
        for i in range(-5, 5 + 1):
            planegrid.add(Line(i * UP + 10 * RIGHT, i * UP + 10 * LEFT))
        planegrid.set_color(GREY)
        planegrid.fade(0.875)
        self.add(planegrid)
        self.wait()
        def tancriccreator(c1, c2, t, p):
            if t %  (1 / 2) == 0:
                t -= 1 / 2048
            h = HomotheticCenter(c1, c2)
            cpt = Line(c2.get_center(), c2.get_start())
            cpt.rotate(t * 2 * PI, about_point = c2.get_center())
            cpt = Line(h[p].get_center(), cpt.get_end())
            cpt = Intersection(cpt, c2).set_color(GREEN)
            ocpt = Line(c2.get_center(), c2.get_start())
            ocpt.rotate(t * 2 * PI, about_point = c2.get_center())
            ocpt = Line(h[p].get_center(), ocpt.get_end())
            ocpt = Intersection(ocpt, c1).set_color(GREEN)
            fcen = Intersection(Line(c2.get_center(), cpt[0].get_center()), Line(c1.get_center(), ocpt[(0 + p) % 2].get_center())).set_color(PINK)
            scen = Intersection(Line(c2.get_center(), cpt[1].get_center()), Line(c1.get_center(), ocpt[(1 + p) % 2].get_center())).set_color(PINK)
            resgrp = VGroup()
            resgrp.add(Circle(arc_center = fcen.get_center(), radius = Line(fcen.get_center(), cpt[0].get_center()).get_length()).set_color(GREEN))
            resgrp.add(Circle(arc_center = scen.get_center(), radius = Line(scen.get_center(), cpt[1].get_center()).get_length()).set_color(GREEN))
            return resgrp
        def circlefamily(c1, c2, p, dt = 0.125):
            resgrp = VGroup()
            for t in np.arange(0, 1, dt):
                resgrp.add(tancriccreator(c1, c2, t, p))
            return resgrp
        c1 = Circle(arc_center = 2.5 * LEFT + DOWN, radius = 0.5)
        c2 = Circle(arc_center = 2.5 * RIGHT + DOWN, radius = 2)
        c3 = Circle(arc_center = 0.5 * LEFT + 1.5 * UP, radius = 1)
        self.play(Write(c1), Write(c2), Write(c3))
        self.wait()
        h12 = HomotheticCenter(c1, c2)
        h23 = HomotheticCenter(c2, c3)
        h31 = HomotheticCenter(c3, c1)
        def extendline(obj):
            return Line(obj.get_start() + 10 * (obj.get_end() - obj.get_start()), obj.get_start() - 10 * (obj.get_end() - obj.get_start()))
        tanlines = VGroup()
        for hcen, cir in zip([h12, h23, h31], [c1, c2, c3]):
            temp = VGroup()
            for obj in hcen:
                tlines = CircleTangentLines(cir, obj)
                for tl in tlines:
                    temp.add(extendline(tl))
            tanlines.add(temp)
        tanlines.set_color(GREEN)
        for obj, hcen in zip(tanlines, [h12, h23, h31]):
            self.play(Write(obj))
            self.play(Write(hcen))
            self.play(FadeOut(obj))
            self.wait()
        homaxes = VGroup(
            extendline(Line(h12[0].get_center(), h23[0].get_center())),
            extendline(Line(h23[0].get_center(), h31[0].get_center())),
            extendline(Line(h31[0].get_center(), h12[0].get_center())),
            extendline(Line(h12[1].get_center(), h23[1].get_center()))
        ).set_color(GOLD)
        self.play(Write(homaxes))
        self.wait()
        homax = TextMobject("Homothetic Axis").move_to(4 * RIGHT + 3 * UP).add_background_rectangle()
        self.play(Write(homax))
        self.wait()
        sols = ApolloniusCircles(c1, c2, c3)
        sols = VGroup(sols[6], sols[7])
        #sols = VGroup(sols[4], sols[5])
        self.play(FadeOut(homax), Write(sols), FadeOut(h12), FadeOut(h23), FadeOut(h31), FadeOut(homaxes))
        self.wait()
        strt = 16
        self.play(c3.fade, 0.5)
        self.wait()
        famcir = circlefamily(c1, c2, 1, dt = 1 / 64)
        for j in range(len(famcir)):
            obj = VGroup(famcir[(j + strt) % 64].fade(0 / 4), famcir[(j - 1 + strt) % 64].fade(1 / 4), famcir[(j - 2 + strt) % 64].fade(2 / 4), famcir[(j - 3 + strt) % 64].fade(3 / 4)).set_color(BLUE)
            self.add(obj)
            self.wait(2 / 15)
            self.remove(obj)
        #self.add(obj)
        self.play(
            Write(h12[1]),
            c3.fade, -1,
        )
        self.wait()
        self.play(c1.fade, 0.5)
        self.wait()
        famcir = circlefamily(c3, c2, 1, dt = 1 / 64)
        for j in range(len(famcir)):
            obj = VGroup(famcir[(j + strt) % 64].fade(0 / 4), famcir[(j - 1 + strt) % 64].fade(1 / 4), famcir[(j - 2 + strt) % 64].fade(2 / 4), famcir[(j - 3 + strt) % 64].fade(3 / 4)).set_color(BLUE)
            self.add(obj)
            self.wait(2 / 15)
            self.remove(obj)
        self.play(
            Write(h23[1]),
            c1.fade, -1,
        )
        self.wait()
        self.play(c2.fade, 0.5)
        self.wait()
        famcir = circlefamily(c3, c1, 1, dt = 1 / 64)
        for j in range(len(famcir)):
            obj = VGroup(famcir[(j + strt) % 64].fade(0 / 4), famcir[(j - 1 + strt) % 64].fade(1 / 4), famcir[(j - 2 + strt) % 64].fade(2 / 4), famcir[(j - 3 + strt) % 64].fade(3 / 4)).set_color(BLUE)
            self.add(obj)
            self.wait(2 / 15)
            self.remove(obj)
        self.play(
            Write(h31[1]),
            c2.fade, -1,
        )
        self.wait()
        homline = extendline(Line(h12[1].get_center(), h23[1].get_center())).set_color(ORANGE)
        self.wait()
        self.play(Write(homline))
        self.wait()
        finstmt = VGroup(
            TextMobject("Any line connecting two Homothetic centers also contains a \\\\ third Homothetic center"),
            TextMobject("The six Homothetic centers are spanned by just four lines \\\\ called the Homothetic Axes"),
            TextMobject("Each Homothetic Axis is a Radical Axis of a pair of circles \\\\ that are mutually tangent to three given circles")
        )
        finstmt.arrange(DOWN, buff = 1)
        finstmt[0][0].set_color(YELLOW)
        finstmt[1][0][3:23].set_color(ORANGE)
        finstmt[1][0][39:48].set_color(BLUE)
        finstmt[1][0][57:].set_color(ORANGE)
        finstmt[2][0][4:18].set_color(ORANGE)
        finstmt[2][0][21:32].set_color(BLUE)
        for txt in finstmt:
            txt.add_background_rectangle()
            self.play(Write(txt))
        self.wait(5)

class FinalSoln(GraphScene, MovingCameraScene):
    CONFIG = {
        "x_min": -10, "x_max": 10, "x_axis_width": 20,
        "y_min": -5, "y_max": 5, "y_axis_height": 10,
        "axes_color": BLACK,
        "graph_origin": ORIGIN,
    }
    def setup(self):
        GraphScene.setup(self)
        MovingCameraScene.setup(self)
    
    def construct(self):
        self.setup_axes()
        self.play(FadeOut(self.axes))
        self.wait()
        planegrid = VGroup()
        for i in range(-10, 10 + 1):
            planegrid.add(Line(5 * UP + i * RIGHT, 5 * DOWN + i * RIGHT))
        for i in range(-5, 5 + 1):
            planegrid.add(Line(i * UP + 10 * RIGHT, i * UP + 10 * LEFT))
        planegrid.set_color(GREY)
        planegrid.fade(0.875)
        self.add(planegrid)
        self.wait()
        c1 = Circle(arc_center = 2.5 * LEFT + DOWN, radius = 0.5)
        c2 = Circle(arc_center = 2.5 * RIGHT + DOWN, radius = 2)
        c3 = Circle(arc_center = 0.5 * LEFT + 1.5 * UP, radius = 1)
        self.play(Write(c1), Write(c2), Write(c3))
        self.wait()
        h12 = HomotheticCenter(c1, c2)
        h23 = HomotheticCenter(c2, c3)
        h31 = HomotheticCenter(c3, c1)
        hcen = VGroup(h12, h23, h31).set_color(GOLD)
        def extendline(obj):
            return Line(obj.get_start() + 10 * (obj.get_end() - obj.get_start()), obj.get_start() - 10 * (obj.get_end() - obj.get_start()))
        homaxes = VGroup(
            extendline(Line(h12[0].get_center(), h23[0].get_center())),
            extendline(Line(h23[0].get_center(), h31[0].get_center())),
            extendline(Line(h31[0].get_center(), h12[0].get_center())),
            extendline(Line(h12[1].get_center(), h23[1].get_center()))
        ).set_color(GOLD)
        r12 = RadicalAxis(c1, c2)
        r23 = RadicalAxis(c2, c3)
        r31 = RadicalAxis(c3, c1)
        raxis = VGroup(r12, r23, r31)
        radcen = RadicalCenter(c1, c2, c3).set_color(MAROON)
        for obj in raxis:
            obj.become(extendline(obj))
            obj.set_color(BLUE)
            self.play(Write(obj))
        self.wait()
        self.play(Write(radcen))
        self.add_foreground_mobject(radcen)
        self.wait()
        self.play(FadeOut(raxis))
        self.wait()
        for obj, cir in zip(hcen, [c1, c2, c3]):
            temp = VGroup()
            temp.add(*CircleTangentLines(cir, obj[0]))
            temp.add(*CircleTangentLines(cir, obj[1]))
            for ob in temp:
                ob.become(extendline(ob))
            temp.set_color(GREEN)
            self.play(Write(temp), Write(obj))
            self.wait()
            self.play(FadeOut(temp))
            self.wait()
        self.play(Write(homaxes))
        self.wait()
        caxis = 2
        fadegrp = VGroup()
        for j in range(4):
            if j == caxis:
                continue
            fadegrp.add(homaxes[j])
        self.play(FadeOut(fadegrp), FadeOut(hcen))
        self.wait()
        discpts = VGroup(
            TextMobject("Homothetic Axis must be a Radical Axis of a \\\\ pair of solution circles"),
            TextMobject("This pair of solution circles are inverses of each other"),
            TextMobject("Points of tangency are also inverses under the same inversion")
        )
        for txt in discpts:
            #txt.scale(0.75)
            txt.add_background_rectangle()
        discpts.arrange(DOWN, buff = 1)
        for txt in discpts:
            self.play(Write(txt))
            self.wait()
        self.play(FadeOut(discpts))
        self.wait(2)
        sols = ApolloniusCircles(c1, c2, c3)
        intpts = VGroup()
        intlines = VGroup()
        tanlines = VGroup()
        tanintersinv = VGroup()
        fadeval = 1 / 2
        for cir in [c1, c2, c3]:
            pt = PerpendicularFoot(Dot(cir.get_center()), homaxes[caxis])
            pt = Circle_Inversion(pt, about_point = cir.get_center(), inversion_radius = Line(cir.get_center(), cir.get_start()).get_length())
            tanintersinv.add(pt)
            pt = Intersection(Line(radcen.get_center(), pt.get_center()), cir)
            intpts.add(pt)
            temp = Line(radcen.get_center(), pt[0].get_center())
            temp = Line(temp.get_start(), temp.get_start() + 10 * (temp.get_end() - temp.get_start()) / temp.get_length()).set_color(BLUE)
            intlines.add(temp)
            temp = VGroup()
            for poin in pt:
                temp.add(extendline(PerpendicularLine(poin, Line(cir.get_center(), poin.get_center()))))
            tanlines.add(temp)
        tanlines.set_color(YELLOW)
        tanlines.fade(fadeval)
        taninters = VGroup()
        for obj in tanlines:
            taninters.add(Intersection(obj[0], obj[1]))
        taninters.fade(fadeval)
        solgrp = VGroup(VGroup(sols[2 * caxis].copy(), sols[2 * caxis + 1].copy()), intpts, intlines)
        for obj in solgrp:
            obj.fade(fadeval)
        self.play(Write(solgrp))
        self.wait()
        self.play(FadeOut(intlines), Write(tanlines[2]), Write(taninters[2]))
        self.wait()
        tanlinesseg = VGroup()
        for i in range(3):
            temp = VGroup()
            temp.add(Line(taninters[i].get_center(), intpts[i][0].get_center()))
            temp.add(Line(taninters[i].get_center(), intpts[i][1].get_center()))
            tanlinesseg.add(temp)
        tanlinesseg.set_color(YELLOW)
        tanlinesseg.fade(fadeval)
        self.play(FadeOut(tanlines[2]), FadeIn(tanlinesseg[2]))
        self.wait()
        self.play(FadeOut(tanlinesseg[2]), FadeOut(taninters[2]))
        self.wait()
        self.play(Write(tanlines[2]))
        self.wait()
        self.play(Write(taninters[2]))
        self.wait()
        self.play(FadeOut(tanlines[2]), FadeIn(tanlinesseg[2]))
        self.wait()
        tanintersinv.set_color(ORANGE)
        #self.play(Write(tanintersinv))
        self.play(
            FadeIn(intlines[2]),
            taninters[2].set_color, BLUE,
        )
        self.wait()
        self.play(Write(tanintersinv[2]),)
        self.wait()
        self.play(
            *[FadeOut(obj) for obj in [solgrp[:2], tanlinesseg[2], intlines[2], taninters[2]]]
        )
        for obj in [intlines, tanlinesseg, taninters, intpts]:
            obj.fade(1 - 1 / (1 - fadeval))
        self.play(Write(tanintersinv[:2]),)
        self.wait()
        self.play(Write(intlines))
        self.wait()
        self.play(Write(intpts))
        self.wait()
        self.play(
            Write(sols[2 * caxis]),
            Write(sols[2 * caxis + 1])
        )
        self.wait()
        self.play(
            *[FadeOut(obj) for obj in [intlines, tanintersinv]]
        )
        self.wait()
        self.play(FadeOut(sols[2 * caxis]), FadeOut(sols[2 * caxis + 1]), FadeOut(intpts), FadeOut(homaxes[caxis]))
        self.wait()
        #starts here...#
        caxis = 0
        intpts = VGroup()
        intlines = VGroup()
        tanlines = VGroup()
        tanintersinv = VGroup()
        for cir in [c1, c2, c3]:
            pt = PerpendicularFoot(Dot(cir.get_center()), homaxes[caxis])
            pt = Circle_Inversion(pt, about_point = cir.get_center(), inversion_radius = Line(cir.get_center(), cir.get_start()).get_length())
            tanintersinv.add(pt)
            pt = Intersection(Line(radcen.get_center(), pt.get_center()), cir)
            intpts.add(pt)
            temp = Line(radcen.get_center(), pt[0].get_center())
            temp = Line(temp.get_start(), temp.get_start() + 10 * (temp.get_end() - temp.get_start()) / temp.get_length()).set_color(BLUE)
            intlines.add(temp)
        tanintersinv.set_color(ORANGE)
        for obj in [homaxes[caxis], tanintersinv, intlines, intpts, VGroup(sols[2 * caxis], sols[2 * caxis + 1])]:
            self.play(Write(obj))
        self.wait()
        self.play(
            *[FadeOut(obj) for obj in [intlines, tanintersinv, sols[2 * caxis], sols[2 * caxis + 1], homaxes[caxis], intpts]]
        )
        self.wait()
        #ends here...#
        #starts here...#
        caxis = 1
        intpts = VGroup()
        intlines = VGroup()
        tanlines = VGroup()
        tanintersinv = VGroup()
        for cir in [c1, c2, c3]:
            pt = PerpendicularFoot(Dot(cir.get_center()), homaxes[caxis])
            pt = Circle_Inversion(pt, about_point = cir.get_center(), inversion_radius = Line(cir.get_center(), cir.get_start()).get_length())
            tanintersinv.add(pt)
            pt = Intersection(Line(radcen.get_center(), pt.get_center()), cir)
            intpts.add(pt)
            temp = Line(radcen.get_center(), pt[0].get_center())
            temp = Line(temp.get_start(), temp.get_start() + 10 * (temp.get_end() - temp.get_start()) / temp.get_length()).set_color(BLUE)
            intlines.add(temp)
        tanintersinv.set_color(ORANGE)
        for obj in [homaxes[caxis], tanintersinv, intlines, intpts, VGroup(sols[2 * caxis], sols[2 * caxis + 1])]:
            self.play(Write(obj))
        self.wait()
        self.play(
            *[FadeOut(obj) for obj in [intlines, tanintersinv, sols[2 * caxis], sols[2 * caxis + 1], homaxes[caxis], intpts]]
        )
        self.wait()
        #ends here...#
        #starts here...#
        caxis = 3
        intpts = VGroup()
        intlines = VGroup()
        tanlines = VGroup()
        tanintersinv = VGroup()
        for cir in [c1, c2, c3]:
            pt = PerpendicularFoot(Dot(cir.get_center()), homaxes[caxis])
            pt = Circle_Inversion(pt, about_point = cir.get_center(), inversion_radius = Line(cir.get_center(), cir.get_start()).get_length())
            tanintersinv.add(pt)
            pt = Intersection(Line(radcen.get_center(), pt.get_center()), cir)
            intpts.add(pt)
            temp = Line(radcen.get_center(), pt[0].get_center())
            temp = Line(temp.get_start(), temp.get_start() + 10 * (temp.get_end() - temp.get_start()) / temp.get_length()).set_color(BLUE)
            intlines.add(temp)
        tanintersinv.set_color(ORANGE)
        for obj in [homaxes[caxis], tanintersinv, intlines, intpts, VGroup(sols[2 * caxis], sols[2 * caxis + 1])]:
            self.play(Write(obj))
        self.wait()
        self.play(
            *[FadeOut(obj) for obj in [intlines, tanintersinv, sols[2 * caxis], sols[2 * caxis + 1], homaxes[caxis], intpts]]
        )
        self.wait()
        #ends here...#
        self.play(
            FadeOut(radcen),
        )
        self.wait()
        for j in range(100):
            self.add(sols[j % 8])
            self.wait(1 / 4)
            self.remove(sols[j % 8])
        self.wait(5)

class GergonneStmt(GraphScene, MovingCameraScene):
    CONFIG = {
        "x_min": -10, "x_max": 10, "x_axis_width": 20,
        "y_min": -5, "y_max": 5, "y_axis_height": 10,
        "axes_color": BLACK,
        "graph_origin": ORIGIN,
    }
    def setup(self):
        GraphScene.setup(self)
        MovingCameraScene.setup(self)
    
    def construct(self):
        self.setup_axes()
        self.play(FadeOut(self.axes))
        self.wait()
        planegrid = VGroup()
        for i in range(-10, 10 + 1):
            planegrid.add(Line(5 * UP + i * RIGHT, 5 * DOWN + i * RIGHT))
        for i in range(-5, 5 + 1):
            planegrid.add(Line(i * UP + 10 * RIGHT, i * UP + 10 * LEFT))
        planegrid.set_color(GREY)
        planegrid.fade(0.875)
        self.add(planegrid)
        self.wait()
        stmt = VGroup(
            TextMobject("Solves the problem with Euclidean construction \\\\ unlike Adrian Van Roomen's"),
            TextMobject("One construction method to solve the whole problem \\\\ unlike Viete's"),
            TextMobject("Minimizes the number of decisions in the constuction \\\\ unlike Newton's")
        ).arrange(DOWN, buff = 1.5)
        stmt[0][0][20:41].set_color(BLUE)
        stmt[0][0][47:].set_color(YELLOW)
        stmt[1][0][:15].set_color(BLUE)
        stmt[1][0][49:].set_color(YELLOW)
        stmt[2][0][12:29].set_color(BLUE)
        stmt[2][0][51:].set_color(YELLOW)
        for txt in stmt:
            txt.add_background_rectangle()
            self.play(Write(txt))
            self.wait()
        self.wait(5)

if __name__ == "__main__":
    module_name = os.path.basename(__file__)
    #command_A = "manim -p  --video_dir ~/Downloads/  "
    clear_cmd = "cls"
    command_A = "python -m manim "
    #command_B = module_name + " " + "HomoCenters" + " -pl -n 90"
    #command_B = module_name + " " + " -pl"
    #command_B = module_name + " " + "ProblemIntroduction" + " -ps -n 16,17"
    #command_B = module_name + " " + "IntroductionScene" + " -p",
    command_B = module_name + " -p"
    #command_B = module_name + " -a"
    #command_B = module_name + " -al"
    os.system(clear_cmd)
    os.system(command_A + command_B)
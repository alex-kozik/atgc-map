#!/usr/bin/python

#################################################################
#                                                               #
#                      MAD MAPPING PROGRAM                      #
#                    PART 3  ( VISUALIZATION )                  #
#                                                               #
#       COPYRIGHT  2004  2005  2006  2007  2008  2009           #
#                        Alexander Kozik                        #
#                                                               #
#                    +----------------------+                   # 
#                    | http://www.atgc.org/ |                   #
#                    +----------------------+                   #
#             +------------------------------------+            #
#             | http://code.google.com/p/atgc-map/ |            #
#             +------------------------------------+            #
#                                                               #
#             UCD Genome Center. R.Michelmore group             #
#                                                               #
#################################################################

#################################################################
#                                   +-------+                   #
#                                   |  BIT  |                   #
#                  SCORING SYSTEM:  |       |                   #
#                                   |  REC  |                   #
#                                   +-------+                   #
#                                                               #
#    .      +-------+-------+-------+-------+-------+-------+   #
#      .    |       |       |       |       |       |       |   #
#        .  |   A   |   B   |   C   |   D   |   H   |   -   |   #
#          .|       |       |       |       |       |       |   #
#   +-------*-------+-------+-------+-------+-------+-------+   #
#   |       | . 6   |  -6   |  -4   |   4   |  -2   |   0   |   #
#   |   A   |       |       |       |       |       |       |   #
#   |       |   0  .|   1   |   1   |   0   |  0.5  |   0   |   #
#   +-------+-------*-------+-------+-------+-------+-------+   #
#   |       |  -6   | . 6   |   4   |  -4   |  -2   |   0   |   #
#   |   B   |       |       |       |       |       |       |   #
#   |       |   1   |   0  .|   0   |   1   |  0.5  |   0   |   #
#   +-------+-------+-------*-------+-------+-------+-------+   #
#   |       |  -4   |   4   | . 4   |  -4   |   0   |   0   |   #
#   |   C   |       |       |       |       |       |       |   #
#   |       |   1   |   0   |   0  .|   1   |   0   |   0   |   #
#   +-------+-------+-------+-------*-------+-------+-------+   #
#   |       |   4   |  -4   |  -4   | . 4   |   0   |   0   |   #
#   |   D   |       |       |       |       |       |       |   #
#   |       |   0   |   1   |   1   |   0  .|   0   |   0   |   #
#   +-------+-------+-------+-------+-------*-------+-------+   #
#   |       |  -2   |  -2   |   0   |   0   | . 2   |   0   |   #
#   |   H   |       |       |       |       |       |       |   #
#   |       |  0.5  |  0.5  |   0   |   0   |   0  .|   0   |   #
#   +-------+-------+-------+-------+-------+-------*-------+   #
#   |       |   0   |   0   |   0   |   0   |   0   | . 0   |   #
#   |   -   |       |       |       |       |       |       |   #
#   |       |   0   |   0   |   0   |   0   |   0   |   0  .|   #
#   +-------+-------+-------+-------+-------+-------+-------*.  #
#                                                               #
#                                                               #
#   NOTES:                                                      #
#      C - NOT A ( H or B )                                     #
#      D - NOT B ( H or A )                                     #
#      H - A and B                                              #
#                                                               #
#################################################################

#################################################################
#                                                               #
#                        EXAMPLES OF SCORING:                   #
#                                                               #
#                                                               #
#  POSITIVE LINKAGE:                                            #
#                                                               #
#  AAAAAAAAAAAAAAAAAAAA     BIT SCORE = 6*20 = 120              #
#  AAAAAAAAAAAAAAAAAAAA     REC SCORE = 0 (0.0)                 #
#                    ..                                         #
#  AAAAAAAAAAAAAAAAAAAA     BIT SCORE = 6*18 - 6*2 = 96         #
#  AAAAAAAAAAAAAAAAAABB     REC SCORE = 2 (2/20 = 0.1)          #
#                                                               #
#  AAAAAAAAAABBBBBBBBBB     BIT SCORE = 6*10 + 6*10 = 120       #
#  AAAAAAAAAABBBBBBBBBB     REC SCORE = 0 (0.0)                 #
#           ..                                                  #
#  AAAAAAAAABABBBBBBBBB     BIT SCORE = 6*18 - 6*2 = 96         #
#  AAAAAAAAAABBBBBBBBBB     REC SCORE = 2 (2/20 = 0.1)          #
#                                                               #
#                                                               #
#  NO LINKAGE:                                                  #
#            ..........                                         #
#  AAAAAAAAAAAAAAAAAAAA     BIT SCORE = 6*10 - 6*10 = 0         #
#  AAAAAAAAAABBBBBBBBBB     REC SCORE = 10 (10/20 = 0.5)        #
#   . . . . . . . . . .                                         #
#  BBBAABBAAAAAAABAABBB     BIT SCORE = 6*10 - 6*10 = 0         #
#  BABBAABBABABABBBAABA     REC SCORE = 10 (10/20 = 0.5)        #
#                                                               #
#                                                               #
#  NEGATIVE LINKAGE:                                            #
#    ..................                                         #
#  AAAAAAAAAAAAAAAAAAAA     BIT SCORE = 6*2 - 6*18 = -96        #
#  AABBBBBBBBBBBBBBBBBB     REC SCORE = 18 (18/20 = 0.9)        #
#    ..................                                         #
#  ABABABABABABABABABAB     BIT SCORE = 6*2 - 6*18 = -96        #
#  ABBABABABABABABABABA     REC SCORE = 18 (18/20 = 0.9)        #
#                                                               #
#                                                               #
#################################################################

def Assign_REC_Color(item_value, rgb_coeff, chick_sat, link_limit):

	rc = 255
	gc = int(round(item_value*rgb_coeff))*chick_sat
	# gc = 255
	if gc >= 255:
		gc = 255
	if gc <= 0:
		gc = 0
	bc = 0
	if item_value >= link_limit:
		rc = 150
		gc = 255
		bc = 150
	if item_value >= 0.49 and item_value <= 0.51:
		rc = 220
		gc = 220
		bc = 220
	if item_value > 0.51:
		rc = 150
		gc = 255
		bc = 255
	if item_value >= 0.60:
		rc = 0
		gc = 200
		bc = 255
	if item_value >= 0.65:
		rc = 0
		gc = 175
		bc = 255
	if item_value >= 0.70:
		rc = 0
		gc = 150
		bc = 255
	if item_value >= 0.75:
		rc = 0
		gc = 100
		bc = 255
	if item_value >= 0.80:
		rc = 0
		gc = 50
		bc = 255
	if item_value >= 0.90:
		rc = 0
		gc = 0
		bc = 255
	#####################
	rvalue = [rc, gc, bc]
	return rvalue

def Assign_BIT_Color(item_value, rgb_coeff, chick_sat, link_limit):

	rc = 255
	# gc = int(round(item_value*rgb_coeff))*chick_sat
	gc = int(round(item_value))
	# gc = 255
	### REVERT GREEN VALUE ###
	gc = 255 - gc + link_limit*2
	if gc >= 255:
		gc = 255
	if gc <= 0:
		gc = 0
	bc = 0
	if item_value <= link_limit:
		rc = 150
		gc = 255
		bc = 150
	if item_value <= 10 and item_value >= -10:
		rc = 220
		gc = 220
		bc = 220
	if item_value < -10:
		rc = 150
		gc = 255
		bc = 255
	if item_value <= -50:
		rc = 0
		gc = 200
		bc = 255
	if item_value <= -100:
		rc = 0
		gc = 150
		bc = 255
	if item_value <= -150:
		rc = 0
		gc = 100
		bc = 255
	if item_value <= -200:
		rc = 0
		gc = 0
		bc = 255
	if item_value <= -300:
		rc = 0
		gc = 0
		bc = 200
	#####################
	bvalue = [rc, gc, bc]
	return bvalue

def Assign_LOD_Color(lvalue):

	print lvalue

def Seqs_Matrix(in_name, list_id, out_name, column_n, diag_fill, dummy_fill, \
		rgb_coeff, chick_sat, link_limit, purple_list, red_list, loc_file, \
		data_type, max_ril_n, link_cut, cell_size, graph_opt, ril_type):

	global cgpdb_style

	print "INPUT FILE:     " + in_name
	print "LIST ID:        " + list_id
	print "LOC FILE:       " + loc_file
	print "OUTPUT FILE:    " + out_name
	print "COLUMN N:       " + `column_n`
	print "DIAG FILL:      " + `diag_fill`
	print "DUMMY FILL:     " + dummy_fill
	print "RGB COEFF:      " + `rgb_coeff`
	print "YELLOW SAT:     " + `chick_sat`
	print "LINKG LIMIT:    " + `link_limit`
	print "LINKAGE CUTOFF: " + `link_cut`
	print "GRAPH OPTION :  " + graph_opt
	print "IMAGE SCALE:    " + cell_size
	print "CROSS TYPE:     " + ril_type

	in_file     = open(in_name,  "rb")
	id_file     = open(list_id,  "rb")
	out_file1   = open(out_name + '.matrix2d.tab', "wb")
	if cgpdb_style == "FALSE":
		out_file2   = out_name + '.matrix2d.large.png'
	### DUMMY CGPDB ###
	if cgpdb_style == "TRUE":
		out_file2   = out_name + '.large.png'

	if graph_opt == "GRAPH":
		out_matrix  = out_name + '.xsparse.matrix'
		out_graph_name = out_name + '.xgraph.large.png'
		out_graph_name2 = out_name + '.xgraph.2500.png'
		out_graph_name3 = out_name + '.xgraph.small.png'

	id_list = []
	loc_id_list = []
	sorted_list = []
	id_array = {}
	loc_id_array = {}
	map_array = {}
	sorted_array = {}
	matrix_array = {}
	frame_array = {}
	high_array = {}
	sequence_array = {}
	ril_array = {}
	ril_list = []
	scores_array = {}
	double_rec = {}
	double_ril = {}
	graph_matrix = {}
	graph_matrix_neg = {}
	node_maxlink = {}
	node_maxlink_neg = {}
	pointX = {}
	pointY = {}
	pointTX = {}
	pointTY = {}
	pointX2 = {}
	pointY2 = {}
	pointTX2 = {}
	pointTY2 = {}

	print_frame = "FALSE"
	### TRY TO READ FRAME MARKERS LIST ###
	try:
		frame_file = open(purple_list)
		print "USING FRAME MARKERS LIST"
		while 1:
			u = frame_file.readline()
			if u == '':
				break
			if '\n' in u:
				u = u[:-1]
			if '\r' in u:
				u = u[:-1]
			u = u.split('\t')
			fm = u[1]
			fl = u[0]
			frame_array[fm] = fl
			print_frame = "TRUE"
	except:
		print "DID NOT FIND FILE:  " + purple_list
		# continue

	print_high = "FALSE"
	try:
		high_file = open(red_list)
		print "USING MARKERS LIST TO HIGHLIGHT"
		while 1:
			g = high_file.readline()
			if g == '':
				break
			if '\n' in g:
				g = g[:-1]
			if '\r' in g:
				g = g[:-1]
			g = g.split('\t')
			hm = g[0]
			high_array[hm] = hm
			print_high = "TRUE"
	except:
		print "DID NOT FIND FILE:  " + red_list

	work_with_loc = "FALSE"
	### TRY TO OPEN LOC FILE #####
	try:
		loc_file_open = open(loc_file)
		work_with_loc = "TRUE"
	except:
		print "DID NOT FIND FILE:  " + loc_file
		draw_ril_image = "FALSE"

	###############################################
	###          PROCESSING LOC FILE            ###

	if work_with_loc == "TRUE":

		### READ LOCUS FILE ###
		print "============================================="
		n = 0
		d = 0
		l = 1
		pool_AB  = {}
		pool_A   = {}
		pool_B   = {}
		markers_A = {}
		markers_B = {}
		ratio_AB = {}
		rils_A   = {}
		rils_B   = {}
		init_len = 10000
		duplic_id = []
		draw_ril_image = "FALSE"

		#################
		while 1:
			t = loc_file_open.readline()
			if t == '':
				break
			if '\n' in t:
				t = t[:-1]
			if '\r' in t:
				t = t[:-1]
			tl = t.split('\t')
			curr_len = len(tl)
			if tl[0][0] == ";":
				print "============================================="
				print tl
				print "                JUNK LINE                    "
				print "============================================="
				time.sleep(2)
			### READ RIL IDs ###
			if len(tl) >= 12:
				if tl[0] == ";" and tl[1] == "1" and tl[2] == "2":
					print "============================================="
					print tl
					print "          RIL IDENTIFIER LINE                "
					print "============================================="
					for ril in tl:
						if ril != ";":
							ril_array[ril] = ril
							ril_list.append(ril)
							ril_int = int(ril)
							if ril_int >= max_ril_n:
								max_ril_n = ril_int
					ril_len = len(tl)
					# draw_ril_image = "TRUE"
					print "MAX RIL NUMBER: " + `max_ril_n`
					time.sleep(2)
					### SET ZERO FOR RILs ####
					for ril_id in ril_list:
						rils_A[ril_id] = 0
						rils_B[ril_id] = 0
					##########################
			if l == 1 and curr_len >= 12 and tl[0][0] != ";":
				init_len = curr_len
			if curr_len == init_len and tl[0][0] != ";":
				####
				id = tl[0]
				# if id not in id_list:
				try:
					id_test = loc_id_array[id]
					duplic_id.append(id)
					print '\n'
					print id + "    IS DUPLICATED -=- CHECK DATA FOR DUPLICATION"
					print ""
					print "REMOVE DUPLICATES FROM LOCUS FILE, THEN CONTINUE"
					print ""
					time.sleep(4)
					d = d + 1
					sys.exit()
				except:
					markers_A[id] = 0
					markers_B[id] = 0
					count_A = 0
					count_B = 0
					loc_id_array[id] = 1
					loc_id_list.append(id)
					q = 1
					while q < init_len:
						data_point = tl[q]
						ril_id = ril_list[q-1]
						if data_point == "A" or data_point == "B" or data_point == "C" \
							or data_point == "D" or data_point == "H" or data_point == "-":
							sequence_array[id,q] = data_point
							scores_array[id,ril_id] = data_point
							if data_point == "A":
								count_A = count_A + 1
								# rils_A[ril_id] = rils_A[ril_id] + 1
							if data_point == "B":
								count_B = count_B + 1
								# rils_B[ril_id] = rils_B[ril_id] + 1
							if data_point == "C" and ril_type == "RIL":
								count_B = count_B + 1
								# rils_B[ril_id] = rils_B[ril_id] + 1
							if data_point == "D" and ril_type == "RIL":
								count_A = count_A + 1
								# rils_A[ril_id] = rils_A[ril_id] + 1
						else:
							print '\n'
							print "WRONG DATA FORMAT"
							print "CHECK LINE:  " + `l` + '\n'
							print t
							print "============================================="
							print "VALUE:   " + data_point
							print "============================================="
							print "    . . .  TERMINATED ..........             "
							print "============================================="
							sys.exit()
						q = q + 1
					sys.stdout.write(".")
					n = n + 1
					count_ALL = count_A + count_B
					markers_A[id] = count_A
					markers_B[id] = count_B
					if count_ALL >= 20:
						pool_AB[id] = int(round(((count_A - count_B)*1.0/count_ALL)*100))
						pool_A[id]  = int(round((count_A*1.0/count_ALL)*100))
						pool_B[id]  = int(round((count_B*1.0/count_ALL)*100))
					if count_ALL < 20:
						pool_AB[id] = "no_data"
				l = l + 1
			if curr_len != init_len and l > 1:
				print '\n'
				print "WRONG NUMBER OF DATA POINTS"
				print "CHECK LINE:  " + `l` + '\n'
				print t
				print "============================================="
				print "   . .     TERMINATED .......... .           "
				print "============================================="
				sys.exit()

		############################
		print ""
		print ril_array
		print ""
		print ril_list
		if len(ril_list)+1 == curr_len:
			draw_ril_image = "TRUE"
			print ""
			print "WILL DRAW RIL IMAGE -=- " + draw_ril_image
			print ""
		if len(ril_list)+1 != curr_len:
			print ""
			print "DO NOT DRAW RIL IMAGE -=- " + draw_ril_image
			print ""
		time.sleep(2)
		############################

		print '\n'
		print "============================================="
		print `n` + " UNIQ IDs IN THE LOCUS FILE FOUND"
		print `d` + " IDs ARE DUPLICATED"

		duplic_id.sort()
		print duplic_id

		print "CONTINUE ANALYSIS WITH " + `n` + " SEQUENCES OUT OF " + `n+d`
		print "============================================="
		print "                SUCCESS!!!                   "
		print "============================================="

		loc_file_open.close()

		###            END OF LOC FILE              ###
		###############################################

	### DRAW OR NOT TO DRAW MAP ###++>
	map_status = "TRUE"
	map_max    = 0
	###############################++<

	font = ImageFont.load_default()

	###  READ SORTED IDs   ###

	while 1:
		t = id_file.readline()
		if t == '':
			break
		if '\n' in t:
			t = t[:-1]
		if '\r' in t:
			t = t[:-1]
		t = t.split('\t')
		test_t = len(t)
		if test_t == 2:
			id = t[0]
		if test_t >= 3:
			id = t[1]
		######  EXTRACT DATA FOR MAP  ######++>
		try:
			if test_t == 2:
				map_array[id] = float(t[1])
			if test_t >= 3:
				map_array[id] = float(t[2])
			if map_array[id] >= map_max:
				map_max = map_array[id]
		except:
			map_array[id] = "NO_DATA"
			map_status = "FALSE"
		####################################++<
		try:
			id_test = sorted_array[id]
			print ""
			print "DUPLICATES: " + id
			print ""
			print "REMOVE DUPLICATES FROM MAP FILE, THEN CONTINUE"
			print ""
			time.sleep(4)
			sys.exit()
		except:
			sorted_array[id] = 1
			sorted_list.append(id)

	### READ ALL HITS FILE ###

	while 1:
		t = in_file.readline()
		if t == '':
			break
		if '\n' in t:
			t = t[:-1]
		if '\r' in t:
			t = t[:-1]
		t = t.split('\t')
		####
		# if t[0] != ";" and t != "" and t[0:4] != 'name':
		if t[0] != ";" and t != "" and t[0:4] != 'name' and len(t) >= 3:
			id_a = t[0]
			id_b = t[1]
			matrix_values = t[column_n + 1]
			# print matrix_values
			# if id_a not in id_list:
			try:
				id_test = id_array[id_a]
			except:
				id_array[id_a] = 1
				id_list.append(id_a)
			try:
				id_test = id_array[id_b]
			except:
				id_array[id_b] = 1
				id_list.append(id_b)

			matrix_array[id_a,id_b] = matrix_values
			matrix_array[id_b,id_a] = matrix_values

	########################################################

	id_list.sort()
	print id_list
	print "------------------------------"
	print "FOUND " + `len(id_list)` + " UNIQ IDs IN MATRIX FILE"
	print "------------------------------"
	print sorted_list
	print "------------------------------"
	print "FOUND " + `len(sorted_list)` + " UNIQ IDs IN MAP FILE"
	print "------------------------------"
	#####################################++>
	print "MAX POSITION:  " + `map_max`
	print "MAP STATUS:    " + map_status
	#####################################++<

	########################################################
	###        MULTI-POINT SPARSE GRAPH MATRIX           ###
	########################################################

	if data_type == "REC" and graph_opt == "GRAPH":
		second_graph = "FALSE"
		###  POSITIVE LINKAGE ###
		for node_A in sorted_list:
			node_maxlink[node_A] = 0.5
			for node_B in sorted_list:
				try:
					linkage = graph_matrix[node_A,node_B]
				except:
					try:
						recombination = float(matrix_array[node_A,node_B])
						linkage = 1 - recombination
						if linkage >= node_maxlink[node_A]:
							node_maxlink[node_A] = linkage
					except:
						recombination = 0.5
						linkage = 1 - recombination
					# if linkage >= 0.6:
					try:
						graph_matrix[node_B,node_A]
					except:
						if linkage >= link_cut:
							graph_matrix[node_A,node_B] = linkage
		###  NEGATIVE LINKAGE ###
		for node_A in sorted_list:
			node_maxlink_neg[node_A] = 0.5
			for node_B in sorted_list:
				try:
					linkage = graph_matrix_neg[node_A,node_B]
				except:
					try:
						recombination = float(matrix_array[node_A,node_B])
						linkage = 1 - recombination
						if linkage <= node_maxlink_neg[node_A]:
							node_maxlink_neg[node_A] = linkage
					except:
						recombination = 0.5
						linkage = 1 - recombination
					# if linkage >= 0.6:
					try:
						graph_matrix_neg[node_B,node_A]
					except:
						if linkage <= 1 - link_cut:
							graph_matrix_neg[node_A,node_B] = linkage

		##################################
		out_file3 = open(out_matrix, "wb")
		### POSITIVE LINKAGE ###
		for node_A in sorted_list:
			for node_B in sorted_list:
				try:
					linkage = graph_matrix[node_A,node_B]
					lvalue = str(round(linkage,2))
					out_file3.write(node_A + '\t' + node_B + '\t' + lvalue + '\n')
				except:
					mooba = "WHATEVER"
		### NEGATIVE LINKAGE ###
		for node_A in sorted_list:
			for node_B in sorted_list:
				try:
					linkage = graph_matrix_neg[node_A,node_B]
					lvalue = str(round(linkage,2))
					out_file3.write(node_A + '\t' + node_B + '\t' + lvalue + '\n')
				except:
					mooba = "WHATEVER"
		# for node in sorted_list:
		#	out_file3.write(node + '\t' + `node_maxlink[node]` + '\n')
		out_file3.close()
		##################################
		### GRAPH IMAGE ###
		graph_k = 20
		graph_margin = 300
		graph_pixels = len(sorted_list)*graph_k
		if graph_pixels <= 1500:
			graph_pixels = 1500
		if graph_pixels >= 10000:
			graph_pixels = 10000
		if graph_pixels > 2500:
			second_graph = "TRUE"
			graph_pixels2 = 2500
			graph_radius2 = 1250
		if graph_pixels >= 10000:
			second_graph = "TRUE"
			graph_pixels2 = 5000
			graph_radius2 = 2500
		graph_radius = graph_pixels/2
		nodes_range = len(sorted_list)
		inner_circle = 260
		mx  = 150
		my =  150
		cs  = 24
		hcs = 12
		txm = 75

		graph_image_size = graph_pixels+graph_margin+cs

		graph_image = Image.new("RGB", ((graph_image_size), (graph_image_size)), (255, 255, 255))
		draw_graph = ImageDraw.Draw(graph_image)

		if second_graph == "TRUE":

			graph_image_size2 = graph_pixels2+graph_margin+cs
			graph_image2 = Image.new("RGB", ((graph_image_size2), (graph_image_size2)), (255, 255, 255))
			draw_graph2 = ImageDraw.Draw(graph_image2)

		tick = 1
		radians_per_degree = 3.141593/180

		#### DRAW NODES (MARKERS)
		for marker in sorted_list:
			tick_update = math.fmod(tick,2)
			angle = tick*(360.0001/nodes_range)
			if tick_update == 0:
				inc = 0
				tx =  txm*1
			if tick_update != 0:
				inc = inner_circle
				tx = -txm*1
			graph_radius = graph_pixels/2 - inc
			graph_r_text = graph_pixels/2 - inc + tx
			########################################
			if second_graph == "TRUE":
				graph_radius2 = graph_pixels2/2 - inc
				graph_r_text2 = graph_pixels2/2 - inc + tx
				############################################
				pointX2[marker] = int(round(graph_radius2*math.sin(angle*radians_per_degree))) + graph_radius2 + inc
				pointY2[marker] = int(round(graph_radius2*math.cos(angle*radians_per_degree))) + graph_radius2 + inc
				pointTX2[marker] = int(round(graph_r_text2*math.sin(angle*radians_per_degree))) + graph_r_text2 + inc - tx
				pointTY2[marker] = int(round(graph_r_text2*math.cos(angle*radians_per_degree))) + graph_r_text2 + inc - tx
				#######################
				draw_graph2.ellipse([pointX2[marker]+mx,pointY2[marker]+my,\
					(pointX2[marker]+mx+cs),(pointY2[marker]+my+cs)],fill=(250,250,0),outline=(0,0,250))
				draw_graph2.text([pointTX2[marker]+mx-hcs,pointTY2[marker]+my],\
					fill=(120,0,200),text=(marker),font=font)
				draw_graph2.line([pointX2[marker]+mx+hcs,pointY2[marker]+my+hcs,\
					(pointTX2[marker]+mx+hcs),(pointTY2[marker]+my+hcs)],fill=(240,0,0))
			####################################################
			pointX[marker] = int(round(graph_radius*math.sin(angle*radians_per_degree))) + graph_radius + inc
			pointY[marker] = int(round(graph_radius*math.cos(angle*radians_per_degree))) + graph_radius + inc
			pointTX[marker] = int(round(graph_r_text*math.sin(angle*radians_per_degree))) + graph_r_text + inc - tx
			pointTY[marker] = int(round(graph_r_text*math.cos(angle*radians_per_degree))) + graph_r_text + inc - tx
			#######################
			draw_graph.ellipse([pointX[marker]+mx,pointY[marker]+my,\
				(pointX[marker]+mx+cs),(pointY[marker]+my+cs)],fill=(250,250,0),outline=(0,0,250))
			draw_graph.text([pointTX[marker]+mx-hcs,pointTY[marker]+my],\
				fill=(120,0,200),text=(marker),font=font)
			draw_graph.line([pointX[marker]+mx+hcs,pointY[marker]+my+hcs,\
				(pointTX[marker]+mx+hcs),(pointTY[marker]+my+hcs)],fill=(240,0,0))

			tick = tick + 1

		#### DRAW EDGES (LINKS)
		tick = 1
		for node_A in sorted_list:
			tick_update = math.fmod(tick,2)
			angle = tick*(360.0001/nodes_range)
			if tick_update == 0:
				inc = 0
			if tick_update != 0:
				inc = inner_circle
			graph_radius = graph_pixels/2 - inc
			if second_graph == "TRUE":
				graph_radius2 = graph_pixels2/2 - inc
			for node_B in sorted_list:
				try:
					linkage = graph_matrix[node_A,node_B]
					if linkage >= 0.9:
						draw_graph.line([(pointX[node_A]+mx+hcs),(pointY[node_A]+my+hcs),\
							(pointX[node_B]+mx+hcs),(pointY[node_B]+my+hcs)],fill=(0,0,0))
						if second_graph == "TRUE":
							draw_graph2.line([(pointX2[node_A]+mx+hcs),(pointY2[node_A]+my+hcs),\
							(pointX2[node_B]+mx+hcs),(pointY2[node_B]+my+hcs)],fill=(0,0,0))
					if linkage < 0.9 and linkage >= 0.85:
						draw_graph.line([(pointX[node_A]+mx+hcs),(pointY[node_A]+my+hcs),\
							(pointX[node_B]+mx+hcs),(pointY[node_B]+my+hcs)],fill=(120,120,120))
						if second_graph == "TRUE":
							draw_graph2.line([(pointX2[node_A]+mx+hcs),(pointY2[node_A]+my+hcs),\
							(pointX2[node_B]+mx+hcs),(pointY2[node_B]+my+hcs)],fill=(120,120,120))
					if linkage < 0.85:
						draw_graph.line([(pointX[node_A]+mx+hcs),(pointY[node_A]+my+hcs),\
							(pointX[node_B]+mx+hcs),(pointY[node_B]+my+hcs)],fill=(180,180,180))
						if second_graph == "TRUE":
							draw_graph2.line([(pointX2[node_A]+mx+hcs),(pointY2[node_A]+my+hcs),\
							(pointX2[node_B]+mx+hcs),(pointY2[node_B]+my+hcs)],fill=(180,180,180))
					if linkage < 0.25:
						draw_graph.line([(pointX[node_A]+mx+hcs),(pointY[node_A]+my+hcs),\
							(pointX[node_B]+mx+hcs),(pointY[node_B]+my+hcs)],fill=(250,0,0))
						if second_graph == "TRUE":
							draw_graph2.line([(pointX2[node_A]+mx+hcs),(pointY2[node_A]+my+hcs),\
							(pointX2[node_B]+mx+hcs),(pointY2[node_B]+my+hcs)],fill=(250,0,0))
				except:
					mooba = "WHATEVER"
			#######################

		if second_graph == "TRUE":

			draw_graph2.text([200,(graph_image_size2-60)],fill=(0,0,0),text=("MAP: " + list_id),font=font)
			draw_graph2.text([600,(graph_image_size2-60)],fill=(0,0,0),text=("MATRIX: " + in_name),font=font)
			linkage_limit2 = round(link_cut,2)
			linkage_limit2 = str(linkage_limit2)
			draw_graph2.text([1000,(graph_image_size2-60)],fill=(0,0,0),text=("LINKAGE CUTOFF: " + linkage_limit2),font=font)

		draw_graph.text([200,(graph_image_size-60)],fill=(0,0,0),text=("MAP: " + list_id),font=font)
		draw_graph.text([600,(graph_image_size-60)],fill=(0,0,0),text=("MATRIX: " + in_name),font=font)
		linkage_limit = round(link_cut,2)
		linkage_limit = str(linkage_limit)
		draw_graph.text([1000,(graph_image_size-60)],fill=(0,0,0),text=("LINKAGE CUTOFF: " + linkage_limit),font=font)

		print "=================================="
		print "      PROCESSING GRAPH IMAGE      "
		print "=================================="
		graph_image.save(out_graph_name, graph_image.format)
		### SMALL IMAGE ####
		mini_graph = Image.open(out_graph_name)
		mini_graph = mini_graph.resize((200,200),Image.ANTIALIAS)
		mini_graph = mini_graph.filter(ImageFilter.SHARPEN)
		mini_graph.save(out_graph_name3, mini_graph.format)
		####################
		if second_graph == "TRUE":
			graph_image2.save(out_graph_name2, graph_image2.format)
		##################################
		# sys.exit()

	#########################################################

	# xn = len(id_list)
	xn = len(sorted_list)

	# size of cell in pixels
	if cell_size == "LARGE":
		c = 10
	if cell_size == "SMALL":
		c = 6
	# margin
	ml = 75		; # LEFT MARGIN
	mr = 200	; # RIGHT MARGIN

	#####################################==>
	# INCREASE RIGHT MARGIN IN THE CASE OF MAP DRAWING
	m_map = 0
	#####################################==<

	#####################################++>
	if map_status == "TRUE":
		map_unit = c*xn/map_max
		print "MAP UNIT:  " + `map_unit` + "  PIXELS"
		### MARGIN FOR MAP
		m_map = 200
		# m_map = 260
		print "MARGIN FOR MAP: " + str(m_map)
	ab_margin = 0
	if work_with_loc == "TRUE":
		ab_margin = 150
	my_lovely_image = Image.new("RGB", ((xn*c+ml+mr+m_map), (xn*c+ml+ml+ab_margin)), (255, 255, 255))
	#####################################==<
	draw_png = ImageDraw.Draw(my_lovely_image)

	######## SET ZERO FOR RIL DATA ########
	for ril_id in ril_list:
		double_ril[ril_id] = 0

	######## RIL IMAGE #########
	if draw_ril_image == "TRUE":
			xrn = max_ril_n
			if len(ril_list) >= max_ril_n:
				xrn = len(ril_list)
			### IMAGE FILES ###
			my_ril_image = Image.new("RGB", ((xrn*c+ml+mr+m_map), (xn*c+ml+ml+ab_margin)), (255, 255, 255))
			draw_ril = ImageDraw.Draw(my_ril_image)
			out_ril_large_name = out_name + '.recomb.large.png'
			out_ril_small_name = out_name + '.recomb.small.png'
			###  TEXT FILES ###
			out_file5 = open(out_name + '.recomb.xmarkers.tab', "wb")	# MARKERS
			out_file6 = open(out_name + '.recomb.xrils.tab', "wb")		# RILs
			out_file5.write("MARKER" + '\t' + "_A_" + '\t' + "_B_" + '\t' + "ALL" + '\t' + "A-B" + '\t' + "---" + '\t' + "DOUBLE" + '\n')
			out_file6.write("RIL_ID" + '\t' + "_A_" + '\t' + "_B_" + '\t' + "ALL" + '\t' + "A-B" + '\t' + "---" + '\t' + "DOUBLE" + '\n')
			### DISPLAY SCORES BY DIFFERENT COLOR ###
			draw_ril.text([ 60,22],fill=(0,0,0),text=("COLOR SCHEME: "),font=font)
			draw_ril.rectangle([155,20,170,35],fill=(250,0,0),outline=(0,0,0))
			draw_ril.text([175,22],fill=(0,0,0),text=("- A"),font=font)
			draw_ril.rectangle([255,20,270,35],fill=(0,0,250),outline=(0,0,0))
			draw_ril.text([275,22],fill=(0,0,0),text=("- B"),font=font)
			draw_ril.rectangle([355,20,370,35],fill=(0,200,0),outline=(0,0,0))
			draw_ril.text([375,22],fill=(0,0,0),text=("- C : not A"),font=font)
			draw_ril.rectangle([455,20,470,35],fill=(250,175,0),outline=(0,0,0))
			draw_ril.text([475,22],fill=(0,0,0),text=("- D : not B"),font=font)
			draw_ril.rectangle([555,20,570,35],fill=(250,250,75),outline=(0,0,0))
			draw_ril.text([575,22],fill=(0,0,0),text=("- H : A+B"),font=font)
			draw_ril.rectangle([655,20,670,35],fill=(120,120,120),outline=(0,0,0))
			draw_ril.text([675,22],fill=(0,0,0),text=("- no data"),font=font)
			#########################################
			y = 0
			for marker in sorted_list:
				# print "XXXXXXXXXXXXX  " + marker + '\t' + `y`
				x = 0
				z = 1
				ro = 0
				go = 0
				bo = 0
				double_rec[marker] = 0
				draw_diag = "FALSE"
				uniform_up = "FALSE"
				uniform_dn = "FALSE"
				while z <= xrn:
				# for ril_id in ril_list:
					ril_id = str(z)
					try:
						score = scores_array[marker,ril_id]
						# double_rec[marker] = 0
						####################################
						### TEST FOR DOUBLE RECOMBINANTS ###
						if y >= 1 and y < len(sorted_list)-1:
							marker_up = sorted_list[y-1]
							marker_dn = sorted_list[y+1]
							ro = 0
							go = 0
							bo = 0
							draw_diag = "FALSE"
							uniform_up = "FALSE"
							uniform_dn = "FALSE"
							try:
								score_up = scores_array[marker_up,ril_id]
								score_dn = scores_array[marker_dn,ril_id]
								#### B CASE ####
								if score_up == "A" and score == "B" and score_dn == "A":
									ro = 250
									go = 250
									bo = 250
									double_rec[marker] = double_rec[marker] + 1
									double_ril[ril_id] = double_ril[ril_id] + 1
									draw_diag = "TRUE"
								if score_up == "D" and score == "B" and score_dn == "D":
									ro = 250
									go = 250
									bo = 250
									double_rec[marker] = double_rec[marker] + 1
									double_ril[ril_id] = double_ril[ril_id] + 1
									draw_diag = "TRUE"
								if score_up == "A" and score == "B" and score_dn == "D":
									ro = 250
									go = 250
									bo = 250
									double_rec[marker] = double_rec[marker] + 1
									double_ril[ril_id] = double_ril[ril_id] + 1
									draw_diag = "TRUE"
								if score_up == "D" and score == "B" and score_dn == "A":
									ro = 250
									go = 250
									bo = 250
									double_rec[marker] = double_rec[marker] + 1
									double_ril[ril_id] = double_ril[ril_id] + 1
									draw_diag = "TRUE"
								#### A CASE ####
								if score_up == "B" and score == "A" and score_dn == "B":
									ro = 250
									go = 250
									bo = 250
									double_rec[marker] = double_rec[marker] + 1
									double_ril[ril_id] = double_ril[ril_id] + 1
									draw_diag = "TRUE"
								if score_up == "C" and score == "A" and score_dn == "C":
									ro = 250
									go = 250
									bo = 250
									double_rec[marker] = double_rec[marker] + 1
									double_ril[ril_id] = double_ril[ril_id] + 1
									draw_diag = "TRUE"
								if score_up == "B" and score == "A" and score_dn == "C":
									ro = 250
									go = 250
									bo = 250
									double_rec[marker] = double_rec[marker] + 1
									double_ril[ril_id] = double_ril[ril_id] + 1
									draw_diag = "TRUE"
								if score_up == "C" and score == "A" and score_dn == "B":
									ro = 250
									go = 250
									bo = 250
									double_rec[marker] = double_rec[marker] + 1
									double_ril[ril_id] = double_ril[ril_id] + 1
									draw_diag = "TRUE"
								#### D CASE ####
								if score_up == "B" and score == "D" and score_dn == "B":
									ro = 250
									go = 250
									bo = 250
									double_rec[marker] = double_rec[marker] + 1
									double_ril[ril_id] = double_ril[ril_id] + 1
									draw_diag = "TRUE"
								if score_up == "C" and score == "D" and score_dn == "C" and ril_type == "RIL":
									ro = 250
									go = 250
									bo = 250
									double_rec[marker] = double_rec[marker] + 1
									double_ril[ril_id] = double_ril[ril_id] + 1
									draw_diag = "TRUE"
								if score_up == "B" and score == "D" and score_dn == "C" and ril_type == "RIL":
									ro = 250
									go = 250
									bo = 250
									double_rec[marker] = double_rec[marker] + 1
									double_ril[ril_id] = double_ril[ril_id] + 1
									draw_diag = "TRUE"
								if score_up == "C" and score == "D" and score_dn == "B" and ril_type == "RIL":
									ro = 250
									go = 250
									bo = 250
									double_rec[marker] = double_rec[marker] + 1
									double_ril[ril_id] = double_ril[ril_id] + 1
									draw_diag = "TRUE"
								#### C CASE ####
								if score_up == "A" and score == "C" and score_dn == "A":
									ro = 250
									go = 250
									bo = 250
									double_rec[marker] = double_rec[marker] + 1
									double_ril[ril_id] = double_ril[ril_id] + 1
									draw_diag = "TRUE"
								if score_up == "D" and score == "C" and score_dn == "D" and ril_type == "RIL":
									ro = 250
									go = 250
									bo = 250
									double_rec[marker] = double_rec[marker] + 1
									double_ril[ril_id] = double_ril[ril_id] + 1
									draw_diag = "TRUE"
								if score_up == "A" and score == "C" and score_dn == "D" and ril_type == "RIL":
									ro = 250
									go = 250
									bo = 250
									double_rec[marker] = double_rec[marker] + 1
									double_ril[ril_id] = double_ril[ril_id] + 1
									draw_diag = "TRUE"
								if score_up == "D" and score == "C" and score_dn == "A" and ril_type == "RIL":
									ro = 250
									go = 250
									bo = 250
									double_rec[marker] = double_rec[marker] + 1
									double_ril[ril_id] = double_ril[ril_id] + 1
									draw_diag = "TRUE"
								#### H CASE ####
								if score_up == "A" and score == "H" and score_dn == "A":
									ro = 120
									go = 120
									bo = 120
									double_rec[marker] = double_rec[marker] + 1
									double_ril[ril_id] = double_ril[ril_id] + 1
									draw_diag = "TRUE"
								if score_up == "B" and score == "H" and score_dn == "B":
									ro = 120
									go = 120
									bo = 120
									double_rec[marker] = double_rec[marker] + 1
									double_ril[ril_id] = double_ril[ril_id] + 1
									draw_diag = "TRUE"
								if score_up == "H" and score == "A" and score_dn == "H":
									ro = 250
									go = 250
									bo = 250
									double_rec[marker] = double_rec[marker] + 1
									double_ril[ril_id] = double_ril[ril_id] + 1
									draw_diag = "TRUE"
								if score_up == "H" and score == "B" and score_dn == "H":
									ro = 250
									go = 250
									bo = 250
									double_rec[marker] = double_rec[marker] + 1
									double_ril[ril_id] = double_ril[ril_id] + 1
									draw_diag = "TRUE"
								#### UNIFORM PATTERN ####
								if score_up == "A" and score == "A":
									uniform_up = "A"
								if score_dn == "A" and score == "A":
									uniform_dn = "A"
								if score_up == "B" and score == "B":
									uniform_up = "B"
								if score_dn == "B" and score == "B":
									uniform_dn = "B"
							############################################################
							except:
								print " -= NO MATCHING PAIR FOUND =- "
						####################################
						if y == len(sorted_list)-1:
							marker_up = sorted_list[y-1]
							ro = 0
							go = 0
							bo = 0
							draw_diag = "FALSE"
							uniform_up = "FALSE"
							uniform_dn = "FALSE"
							try:
								score_up = scores_array[marker_up,ril_id]
								if score_up == "A" and score == "A":
									uniform_up = "A"
								if score_up == "B" and score == "B":
									uniform_up = "B"
							except:
								print " -= NO MATCHING PAIR FOUND =- "
						####################################
						if score == "A":
							colour = [250,0,0]
							draw_ril.rectangle([(ml+x*c),(ml+y*c),(ml+x*c+c),(ml+y*c+c)],fill=(250,0,0),outline=(ro,go,bo))
							rils_A[ril_id] = rils_A[ril_id] + 1
						if score == "B":
							colour = [0,0,250]
							draw_ril.rectangle([(ml+x*c),(ml+y*c),(ml+x*c+c),(ml+y*c+c)],fill=(0,0,250),outline=(ro,go,bo))
							rils_B[ril_id] = rils_B[ril_id] + 1
						if score == "C":
							colour = [0,200,0]	# NOT A
							draw_ril.rectangle([(ml+x*c),(ml+y*c),(ml+x*c+c),(ml+y*c+c)],fill=(0,200,0),outline=(ro,go,bo))
							if ril_type == "RIL":
								rils_B[ril_id] = rils_B[ril_id] + 1
						if score == "D":
							colour = [250,175,0]	# NOT B
							draw_ril.rectangle([(ml+x*c),(ml+y*c),(ml+x*c+c),(ml+y*c+c)],fill=(250,175,0),outline=(ro,go,bo))
							if ril_type == "RIL":
								rils_A[ril_id] = rils_A[ril_id] + 1
						if score == "H":
							colour = [250,250,75]	# A + B
							draw_ril.rectangle([(ml+x*c),(ml+y*c),(ml+x*c+c),(ml+y*c+c)],fill=(250,250,75),outline=(ro,go,bo))
						if score == "-":
							colour = [120,120,120]
							draw_ril.rectangle([(ml+x*c),(ml+y*c),(ml+x*c+c),(ml+y*c+c)],fill=(120,120,120),outline=(ro,go,bo))
						# print score + '\t' + `x` + '\t' + `y`
						if uniform_up == "A":
							draw_ril.line([(ml+x*c+1),(ml+y*c),(ml+x*c+c-1),(ml+y*c)],fill=(200,0,0))
						if uniform_up == "B":
							draw_ril.line([(ml+x*c+1),(ml+y*c),(ml+x*c+c-1),(ml+y*c)],fill=(0,0,200))
						if uniform_dn == "A":
							if y > len(sorted_list)-1:
								draw_ril.line([(ml+x*c+1),(ml+y*c+c),(ml+x*c+c-1),(ml+y*c+c)],fill=(200,0,0))
						if uniform_dn == "B":
							if y > len(sorted_list)-1:
								draw_ril.line([(ml+x*c+1),(ml+y*c+c),(ml+x*c+c-1),(ml+y*c+c)],fill=(0,0,200))
						if draw_diag == "TRUE":
							draw_ril.line([(ml+x*c),(ml+y*c),(ml+x*c+c),(ml+y*c+c)],fill=(ro,go,bo))
							draw_ril.line([(ml+x*c+c),(ml+y*c),(ml+x*c),(ml+y*c+c)],fill=(ro,go,bo))
					except:
						# print "NO DATA FOR " + marker
						draw_ril.rectangle([(ml+x*c),(ml+y*c),(ml+x*c+c),(ml+y*c+c)],fill=(60,60,60),outline=(0,0,0))
					x = x + 1
					z = z + 1
				### MARKER TEXT FILE
				try:
					ab_value_marker = int(round(((markers_A[marker]-markers_B[marker])*1.0/(markers_A[marker]+markers_B[marker]))*100))
					a_value = markers_A[marker]
					b_value = markers_B[marker]
					ab_value = a_value + b_value
					double_marker = double_rec[marker]
				except:
					ab_value_marker = "-"
					a_value = "-"
					b_value = "-"
					ab_value = "-"
					double_marker = "-"
				ab_value_marker = str(ab_value_marker)
				a_value = str(a_value)
				b_value = str(b_value)
				ab_value = str(ab_value)
				double_marker = str(double_marker)
				out_file5.write(marker + '\t' + a_value + '\t' + b_value + '\t' \
				+ ab_value + '\t' + ab_value_marker + '\t' \
				+ "---" + '\t' + double_marker + '\n')
				####################
				y = y + 1
			#########################################
	############################

	q = 0
	################## BIT SCORES #######################
	if data_type == "BIT":
		if cell_size == "LARGE":
			cell_k = 3.0
		if cell_size == "SMALL":
			cell_k = 5.0
		test_list = [300, 250, 200, 150, 100, 50, 10, 0, -10, -50, -100, -150, -200, -250, -300]
		draw_png.text([100,12],fill=(0,0,0),text=("BIT SCORE SCALE: "),font=font)
		for item in test_list:
			# print item
			item_value = float(item)
			#####################################==>
			bvalue = Assign_BIT_Color(item_value-1, rgb_coeff, chick_sat, link_limit)
			rc = bvalue[0]
			gc = bvalue[1]
			bc = bvalue[2]
			#####################################==>
			if item_value != -300:
				draw_png.rectangle([200+q*cell_k*c,15,200+cell_k*c+q*cell_k*c,25],fill=(rc,gc,bc),outline=(0,0,0))
			draw_png.text([190+q*cell_k*c,28],fill=(0,0,124),text=(str(item)),font=font)

			if item_value == -300:
				draw_png.rectangle([200+q*cell_k*c+30,15,200+cell_k*c+q*cell_k*c+30,25],fill=(82,82,82),outline=(0,0,0))
				draw_png.text([190+q*cell_k*c+40,28],fill=(82,82,82),text=("no data"),font=font)
			#####################################==<
			q = q+1

		####### END OF BIT SCORES #######

	################## REC SCORES #######################
	if data_type == "REC":
		if cell_size == "LARGE":
			cell_k = 3.0
		if cell_size == "SMALL":
			cell_k = 5.0
		test_list = [0.00, 0.10, 0.20, 0.30, 0.40, 0.49, 0.50, 0.51, 0.60, 0.65, 0.70, 0.75, 0.80, 0.90, 1.00]
		# draw_png.text([100,12],fill=(0,0,0),text=("DISTANCE SCALE: "),font=font)
		draw_png.text([100,12],fill=(0,0,0),text=("LINKAGE SCALE: "),font=font)
		for item in test_list:
			# print item
			item_value = float(item)
			#####################################==>
			rvalue = Assign_REC_Color(item_value+0.0001, rgb_coeff, chick_sat, link_limit)
			rc = rvalue[0]
			gc = rvalue[1]
			bc = rvalue[2]
			####################
			if item_value != 1.0:
				draw_png.rectangle([(200+q*cell_k*c),15,(200+cell_k*c+q*cell_k*c),25],fill=(rc,gc,bc),outline=(0,0,0))
			draw_png.text([190+q*cell_k*c,28],fill=(0,0,124),text=(str(item)),font=font)

			if item_value == 1.0:
				draw_png.rectangle([(200+q*cell_k*c+30),15,(200+cell_k*c+q*cell_k*c+30),25],fill=(82,82,82),outline=(0,0,0))
				draw_png.text([190+q*cell_k*c+40,28],fill=(82,82,82),text=("no data"),font=font)
			#####################################==<
			q = q+1

		####### END OF REC SCORES #######

	################## LOD SCORES #######################

	if data_type == "LOD":
		if cell_size == "LARGE":
			cell_k = 3.0
		if cell_size == "SMALL":
			cell_k = 5.0
		test_list = [10, 9, 8, 7, 6, 5, 4, 3, 2, 1, 0]
		draw_png.text([100,12],fill=(0,0,0),text=("LOD SCORE SCALE: "),font=font)
		for item in test_list:
			# print item
			item_value = float(item)
			rc = 255
			# gc = int(round(item_value*rgb_coeff))*chick_sat
			gc = int(round(item_value*25))
			### REVERT GREEN VALUE ###
			gc = 255 - gc + link_limit*25
			# gc = 255
			if gc >= 255:
				gc = 255
			if gc <= 0:
				gc = 0
			bc = 0
			if item_value <= link_limit:
				rc = 0
				gc = 255
				bc = 0
			# if item_value >= 0.7:
			if item_value <= 2:
				rc = 0
				gc = 255
				bc = 255
			if item_value <= 1:
				rc = 0
				gc = 164
				bc = 255
			if item_value <= 0:
				rc = 0
				gc = 82
				bc = 255

			#####################################==>
			if item_value != 0:
				draw_png.rectangle([200+q*cell_k*c,15,200+cell_k*c+q*cell_k*c,25],fill=(rc,gc,bc),outline=(0,0,0))
			draw_png.text([190+q*cell_k*c,28],fill=(0,0,124),text=(str(item)),font=font)

			if item_value == 0:
				draw_png.rectangle([200+q*cell_k*c+30,15,200+cell_k*c+q*cell_k*c+30,25],fill=(82,82,82),outline=(0,0,0))
				draw_png.text([190+q*cell_k*c+40,28],fill=(82,82,82),text=("no data"),font=font)
			#####################################==<
			q = q+1

		####### END OF LOD SCORES #######

	k = 0
	out_file1.write("#" + '\t')
	#####################################++>
	if map_status == "TRUE":
		draw_png.text([(mr+xn*c+m_map-270),(ml-c/2-25)],text=('MARKERS:'),font=font,fill=(0,0,124))
		draw_png.text([(mr+xn*c+m_map-25),(ml-c/2-25)],text=('MAP:'),font=font,fill=(0,0,124))
		draw_png.line([(mr+xn*c+m_map-20),(ml-c/2),(mr+xn*c+m_map-20),(ml+map_max*map_unit+c/2)],fill=(0,0,64))
		draw_png.text([(mr+xn*c+m_map-5),(ml-c/2-5)],text=('0'),font=font,fill=(0,0,124))
		draw_png.text([(mr+xn*c+m_map-8),(ml+map_max*map_unit+c/2)],text=(str(round(map_max,1))),font=font,fill=(0,0,124))
		#### DRAW TICKS ON MAP ####
		tic = 0
		while tic <= map_max:
			draw_png.line([(mr+xn*c+m_map-18),(ml+map_unit*tic),(mr+xn*c+m_map-15),(ml+map_unit*tic)],fill=(0,0,64))
			if tic != 0:
				draw_png.text([(mr+xn*c+m_map-8),(ml+map_unit*tic-5)],text=(`tic`),font=font,fill=(0,0,124))
			tic = tic + 10
	######       RIL MAP    ################
	if draw_ril_image == "TRUE":
		draw_ril.text([(mr+xrn*c+m_map-270),(ml-c/2-25)],text=('MARKERS:'),font=font,fill=(0,0,124))
		draw_ril.text([(mr+xrn*c+m_map-25),(ml-c/2-25)],text=('MAP:'),font=font,fill=(0,0,124))
		draw_ril.line([(mr+xrn*c+m_map-20),(ml-c/2),(mr+xrn*c+m_map-20),(ml+map_max*map_unit+c/2)],fill=(0,0,64))
		draw_ril.text([(mr+xrn*c+m_map-5),(ml-c/2-5)],text=('0'),font=font,fill=(0,0,124))
		draw_ril.text([(mr+xrn*c+m_map-8),(ml+map_max*map_unit+c/2)],text=(str(round(map_max,1))),font=font,fill=(0,0,124))
		#### DRAW TICKS ON MAP ####
		tic = 0
		while tic <= map_max:
			draw_ril.line([(mr+xrn*c+m_map-18),(ml+map_unit*tic),(mr+xrn*c+m_map-15),(ml+map_unit*tic)],fill=(0,0,64))
			if tic != 0:
				draw_ril.text([(mr+xrn*c+m_map-8),(ml+map_unit*tic-5)],text=(`tic`),font=font,fill=(0,0,124))
			tic = tic + 10
	### TICKS ON X AXIS OF RIL IMAGE ####
	if draw_ril_image == "TRUE":
		s = 0
		if s == 0:
			draw_ril.text([(ml+s*c),(ml-18)],fill=(0,0,0),text=("RILs:"),font=font)
		while s < xrn:
			tick_update = math.fmod(s+1,5)
			if tick_update == 0:
				move_me_x = 0
				# if s+1 >= 100:
				#	move_me_x == 5
				draw_ril.text([(ml+s*c-move_me_x),(ml-18)],fill=(0,0,0),text=(`s+1`),font=font)
				draw_ril.line([ml+s*c+c/2,ml-2,ml+s*c+c/2,ml-6],fill=(0,0,0))
			s = s + 1
	#####################################++<
	for id in sorted_list:
		#####################################==>
		## TICK ON X AXIS
		tick_update = math.fmod(k+1,5)
		#####################################==<
		if tick_update == 0:
			#####################################==>
			move_me_x = 0
			# if k+1 >= 100:
			#	move_me_x == 5
			draw_png.text([(ml+k*c-move_me_x),(ml-18)],fill=(0,0,0),text=(`k+1`),font=font)
			draw_png.line([ml+k*c+c/2,ml-2,ml+k*c+c/2,ml-6],fill=(0,0,0))
			#####################################==<
		out_file1.write(sorted_list[k])
		if k < xn - 1:
			out_file1.write('\t')
		if k == xn - 1:
			out_file1.write('\n')
		#####################################++>
		### FINE MAP LAYOUT ###
		if map_status == "TRUE":
			tick_shift = math.fmod(k+1,2)
			if cell_size == "SMALL":
				id_shift = 0
				line_shift = 65
				# line_shift = 125
			if cell_size == "LARGE":
				id_shift = 0
				line_shift = 10
			if tick_shift == 0 and cell_size == "SMALL":
				# id_shift = 65
				id_shift = 125
				# line_shift = 0
				line_shift = -60
			if tick_shift == 0 and cell_size == "LARGE":
				id_shift = 0
				line_shift = 0

			### LINES ON MAP ###
			if tick_shift != 0:
				# draw_png.line([(ml+xn*c+m_map-30-line_shift),(ml+k*c+c/2),(ml+xn*c+m_map-10),(ml+k*c+c/2)],fill=(220,120,160))
				draw_png.line([(ml+xn*c+m_map+30-line_shift),(ml+k*c+c/2),(ml+xn*c+m_map+20),(ml+k*c+c/2)],fill=(120,120,164))
				# draw_png.line([(ml+xn*c+m_map-10),(ml+k*c+c/2),(ml+xn*c+m_map+80),(ml+map_array[id]*map_unit)],fill=(120,220,160))
				draw_png.line([(ml+xn*c+m_map+20),(ml+k*c+c/2),(ml+xn*c+m_map+85),(ml+map_array[id]*map_unit)],fill=(120,120,164))
			draw_png.line([(ml+xn*c+m_map+85),(ml+map_array[id]*map_unit),(ml+xn*c+m_map+100),(ml+map_array[id]*map_unit)],fill=(120,120,160))
			#################################

			draw_png.line([(ml+xn*c+m_map-195),(ml+k*c+c/2),(ml+xn*c+m_map-170+id_shift),(ml+k*c+c/2)],fill=(120,120,160))
			draw_png.text([(ml+xn*c+m_map-160+id_shift),(ml+k*c+c/2-5)],fill=(0,0,0),text=(id),font=font)
			if print_frame == "TRUE":
				try:
					frame_id = frame_array[id]
					draw_png.text([(ml+xn*c+m_map-160+id_shift),(ml+k*c+c/2-5)],fill=(0,0,225),text=(id),font=font)
					# if cell_size == "SMALL":
					#	draw_png.text([(ml+xn*c+m_map-160+id_shift/10),(ml+k*c+c/2-5)],fill=(0,0,225),text=(frame_id),font=font)
					#	# draw_png.text([(ml+xn*c+m_map-185+id_shift/10),(ml+k*c+c/2-5)],fill=(0,0,225),text=(frame_id),font=font)
					if cell_size == "LARGE":
						draw_png.text([(ml+xn*c+m_map-15+id_shift),(ml+k*c+c/2-5)],fill=(0,0,225),text=(frame_id),font=font)
						# draw_png.text([(ml+xn*c+m_map-35+id_shift),(ml+k*c+c/2-5)],fill=(0,0,225),text=(frame_id),font=font)
					print id + '\t' + " ++++ IT IS A FRAME MARKER ++++ "
				except:
					# continue
					print id + '\t' + "--- IT IS NOT A FRAME MARKER ---"
			if print_high == "TRUE":
				try:
					high_id = high_array[id]
					draw_png.text([(ml+xn*c+m_map-160+id_shift),(ml+k*c+c/2-5)],fill=(255,0,0),text=(id),font=font)
					# if cell_size == "SMALL":
					#	draw_png.text([(ml+xn*c+m_map-180+id_shift/10),(ml+k*c+c/2-5)],fill=(255,0,0),text=("*"),font=font)
					if cell_size == "LARGE":
						draw_png.text([(ml+xn*c+m_map-25+id_shift),(ml+k*c+c/2-5)],fill=(255,0,0),text=("*"),font=font)
					print id + '\t' + "IT IS A RED MARKER"
				except:
					print id + '\t' + "IT IS NOT A NEW MARKER"
			# draw_png.line([(ml+xn*c+m_map-30-line_shift),(ml+k*c+c/2),(ml+xn*c+m_map-10),(ml+k*c+c/2)],fill=(120,120,160))
			# draw_png.line([(ml+xn*c+m_map-10),(ml+k*c+c/2),(ml+xn*c+m_map+80),(ml+map_array[id]*map_unit)],fill=(120,120,160))
			# draw_png.line([(ml+xn*c+m_map+80),(ml+map_array[id]*map_unit),(ml+xn*c+m_map+100),(ml+map_array[id]*map_unit)],fill=(120,120,160))
		##################################### RIL IMAGE ###############################################
			if draw_ril_image == "TRUE":
				draw_ril.line([(ml+xrn*c+m_map-195),(ml+k*c+c/2),(ml+xrn*c+m_map-170+id_shift),(ml+k*c+c/2)],fill=(120,120,164))
				if tick_shift != 0:
					# draw_ril.text([(ml+xrn*c+m_map-160+id_shift),(ml+k*c+c/2-5)],fill=(0,0,0),text=(id),font=font)
					# draw_ril.line([(ml+xrn*c+m_map-30-line_shift),(ml+k*c+c/2),(ml+xrn*c+m_map-10),(ml+k*c+c/2)],fill=(120,120,164))
					# draw_ril.line([(ml+xrn*c+m_map-10),(ml+k*c+c/2),(ml+xrn*c+m_map+80),(ml+map_array[id]*map_unit)],fill=(120,120,164))
					draw_ril.line([(ml+xrn*c+m_map+30-line_shift),(ml+k*c+c/2),(ml+xrn*c+m_map+20),(ml+k*c+c/2)],fill=(120,120,164))
					draw_ril.line([(ml+xrn*c+m_map+20),(ml+k*c+c/2),(ml+xrn*c+m_map+85),(ml+map_array[id]*map_unit)],fill=(120,120,164))
				draw_ril.line([(ml+xrn*c+m_map+85),(ml+map_array[id]*map_unit),(ml+xrn*c+m_map+100),(ml+map_array[id]*map_unit)],fill=(120,120,164))
				draw_ril.text([(ml+xrn*c+m_map-160+id_shift),(ml+k*c+c/2-5)],fill=(0,0,0),text=(id),font=font)
				#### DOUBLE RECOMBINANATS ####
				trouble = double_rec[id]
				if trouble > 0 and trouble < 2:
					draw_ril.text([(ml+xrn*c+m_map-160+id_shift),(ml+k*c+c/2-5)],fill=(0,60,180),text=(id),font=font)
					# if cell_size == "SMALL":
                                        #        draw_ril.text([(ml+xrn*c+m_map-185+id_shift/10),(ml+k*c+c/2-5)],fill=(0,60,180),text=(`trouble`),font=font)
					if cell_size == "LARGE":
						draw_ril.text([(ml+xrn*c+m_map-15+id_shift),(ml+k*c+c/2-5)],fill=(0,60,180),text=(`trouble`),font=font)
				if trouble == 2:
					draw_ril.text([(ml+xrn*c+m_map-160+id_shift),(ml+k*c+c/2-5)],fill=(0,60,180),text=(id),font=font)
					# if cell_size == "SMALL":
                                        #        draw_ril.text([(ml+xrn*c+m_map-185+id_shift/10),(ml+k*c+c/2-5)],fill=(125,0,0),text=(`trouble`),font=font)
					if cell_size == "LARGE":
						draw_ril.text([(ml+xrn*c+m_map-15+id_shift),(ml+k*c+c/2-5)],fill=(125,0,0),text=(`trouble`),font=font)
				if trouble > 2:
					draw_ril.text([(ml+xrn*c+m_map-160+id_shift),(ml+k*c+c/2-5)],fill=(0,60,180),text=(id),font=font)
					# if cell_size == "SMALL":
					#	draw_ril.text([(ml+xrn*c+m_map-185+id_shift/10),(ml+k*c+c/2-5)],fill=(200,0,0),text=(`trouble`),font=font)
					if cell_size == "LARGE":
						draw_ril.text([(ml+xrn*c+m_map-15+id_shift),(ml+k*c+c/2-5)],fill=(200,0,0),text=(`trouble`),font=font)
				##############################
		###############################################################################################
		#####################################++<
		k = k + 1

	## WRITE IDs ON MAP AGAIN ##
	kkk = 0
	for id in sorted_list:

		if map_status == "TRUE" and cell_size == "SMALL":
                        tick_shift = math.fmod(kkk+1,2)
                        id_shift = 0
                        line_shift = 65
                        # line_shift = 125
                        if tick_shift == 0:
                                # id_shift = 65
                                id_shift = 125
                                # line_shift = 0
                                line_shift = -60

			draw_png.text([(ml+xn*c+m_map-160+id_shift),(ml+kkk*c+c/2-5)],fill=(0,0,0),text=(id),font=font)

			if draw_ril_image == "TRUE":

				draw_ril.text([(ml+xrn*c+m_map-160+id_shift),(ml+kkk*c+c/2-5)],fill=(0,0,0),text=(id),font=font)

		kkk = kkk + 1
	## END OF DUMMY IDs ##

	# y position on image
	yn = 0
	for id in sorted_list:
		k  = 0
		out_file1.write(id + '\t')
		#####################################==>
		tick_update = math.fmod(yn+1,5)
		#####################################==<
		if tick_update == 0:
			move_me_x = 0
			if yn+1 >= 100:
				move_me_x = 5
			if yn+1 >= 1000:
				move_me_x = 10
			#####################################==>
			draw_png.text([(ml-25-move_me_x),(ml+yn*c-3)],fill=(0,0,0),text=(`yn+1`),font=font)
			#####################################==<
			#####################################==>
			draw_png.line([ml-5,ml+yn*c+c/2,ml-2,ml+yn*c+c/2],fill=(0,0,0))
			### Y AXIS ON RIL IMAGE ###
			if draw_ril_image == "TRUE":
				draw_ril.text([(ml-25-move_me_x),(ml+yn*c-3)],fill=(0,0,0),text=(`yn+1`),font=font)
				draw_ril.line([ml-5,ml+yn*c+c/2,ml-2,ml+yn*c+c/2],fill=(0,0,0))
			#####################################==<
		######## BIT SCORE CONDITION ###########
		while k < xn and data_type == "BIT":
			idx = sorted_list[k]
			try:
				current_value = matrix_array[id,idx]
				image_value = current_value
			except:
				current_value = dummy_fill
				# image_value = 0
				# Maria settings
				image_value = 1
			# if id == idx:
			#	current_value = diag_fill
			#	image_value = diag_fill
			# text file
			out_file1.write(str(current_value))
			if k < xn - 1:
				out_file1.write('\t')
			if k == xn - 1:
				out_file1.write('\n')
			# image processing
			# red green blue
			image_value = float(image_value)

			bvalue = Assign_BIT_Color(image_value, rgb_coeff, chick_sat, link_limit)
			rc = bvalue[0]
			gc = bvalue[1]
			bc = bvalue[2]

			# dummy condition
			if current_value == dummy_fill:
				rc = 82
				gc = 82
				bc = 82

			# draw cell
			draw_png.rectangle([(ml+k*c),(ml+yn*c),(ml+k*c+c),(ml+yn*c+c)],fill=(rc,gc,bc),outline=(0,0,0))
			k = k + 1
		####### END OF BIT SCORE CONDITION ########

		####### REC SCORE CONDITION ########
		while k < xn and data_type == "REC":
			idx = sorted_list[k]
			try:
				current_value = matrix_array[id,idx]
				image_value = current_value
			except:
				current_value = dummy_fill
				# image_value = 0
				# Maria settings
				image_value = 1
			if id == idx:
				current_value = diag_fill
				image_value = diag_fill
			# text file
			out_file1.write(str(current_value))
			if k < xn - 1:
				out_file1.write('\t')
			if k == xn - 1:
				out_file1.write('\n')
			# image processing
			# red green blue
			image_value = float(image_value)

			rvalue = Assign_REC_Color(image_value, rgb_coeff, chick_sat, link_limit)
			rc = rvalue[0]
			gc = rvalue[1]
			bc = rvalue[2]

			# dummy condition
			if current_value == dummy_fill:
				rc = 82
				gc = 82
				bc = 82

			# draw cell
			draw_png.rectangle([(ml+k*c),(ml+yn*c),(ml+k*c+c),(ml+yn*c+c)],fill=(rc,gc,bc),outline=(0,0,0))
			k = k + 1
		####### END OF REC SCORE CONDITION ########

		####### LOD SCORE CONDITION ########
		while k < xn and data_type == "LOD":
			idx = sorted_list[k]
			try:
				current_value = matrix_array[id,idx]
				image_value = current_value
			except:
				current_value = dummy_fill
				# image_value = 0
				# Maria settings
				image_value = 1
			# if id == idx:
			#	current_value = diag_fill
			#	image_value = diag_fill
			# text file
			out_file1.write(str(current_value))
			if k < xn - 1:
				out_file1.write('\t')
			if k == xn - 1:
				out_file1.write('\n')
			# image processing
			# red green blue
			image_value = float(image_value)
			rc = 255
			# JM settings
			gc = int(round(image_value*25))
			### REVERT GREEN VALUE ###
			gc = 255 - gc + link_limit*25
			# gc = 255
			if gc >= 255:
				gc = 255
			if gc <= 0:
				gc = 0
			bc = 0
			if image_value <= link_limit:
				rc = 0
				gc = 255
				bc = 0
			# if image_value >= 0.7:
			if image_value <= 2:
				rc = 0
				gc = 255
				bc = 255
			# if image_value >= 0.8:
			if image_value <= 1:
				rc = 0
				gc = 164
				bc = 255
			# if image_value >= 0.9:
			if image_value <= 0:
				rc = 0
				gc = 82
				bc = 255
			# dummy condition
			if current_value == dummy_fill:
				rc = 82
				gc = 82
				bc = 82

			# draw cell
			draw_png.rectangle([(ml+k*c),(ml+yn*c),(ml+k*c+c),(ml+yn*c+c)],fill=(rc,gc,bc),outline=(0,0,0))
			k = k + 1
		####### END OF LOD SCORE CONDITION ########

		yn = yn + 1

	if work_with_loc == "TRUE":
		k = 0
		for id in sorted_list:
			try:
				ab_value = pool_AB[id]
				a_value  = pool_A[id]
				b_value  = pool_B[id]
				ab_diff  = a_value - b_value
			except:
				ab_value = "no_data"
			if ab_value != "no_data":
				if ab_value >=0:
					rc = 255
					gc = 255
					bc = 0
				if ab_value < 0:
					rc = 0
					gc = 255
					bc = 255
				if a_value != 0 and b_value != 0:
					draw_png.rectangle([(ml+k*c),(ml+yn*c+ab_margin/2+10),(ml+k*c+c),\
						(ml+yn*c+ab_margin/2-a_value+10)],fill=(255,0,0),outline=(255,0,0))
					draw_png.rectangle([(ml+k*c),(ml+yn*c+ab_margin/2+10),(ml+k*c+c),\
						(ml+yn*c+ab_margin/2+b_value+10)],fill=(0,0,255),outline=(0,0,255))
					draw_png.rectangle([(ml+k*c),(ml+yn*c+ab_margin/2+10),(ml+k*c+c),\
						(ml+yn*c+ab_margin/2+10-ab_value)],fill=(rc,gc,bc),outline=(rc,gc,bc))
					k = k + 1
				if a_value == 0:
					draw_png.rectangle([(ml+k*c),(ml+yn*c+ab_margin/2),(ml+k*c+c),\
						(ml+yn*c+ab_margin/2+b_value-20)],fill=(0,0,255),outline=(0,0,255))
					draw_png.rectangle([(ml+k*c),(ml+yn*c+ab_margin/2+b_value-17),(ml+k*c+c),\
						(ml+yn*c+ab_margin/2+b_value-12)],fill=(0,0,255),outline=(0,0,255))
					k = k + 1
				if b_value == 0:
					draw_png.rectangle([(ml+k*c),(ml+yn*c+ab_margin/2+20),(ml+k*c+c),\
						(ml+yn*c+ab_margin/2-a_value+40)],fill=(255,0,0),outline=(255,0,0))
					draw_png.rectangle([(ml+k*c),(ml+yn*c+ab_margin/2-a_value+37),(ml+k*c+c),\
						(ml+yn*c+ab_margin/2-a_value+32)],fill=(255,0,0),outline=(255,0,0))

					k = k + 1
			if ab_value == "no_data":
				rc = 82
				gc = 82
				bc = 82
				draw_png.rectangle([(ml+k*c),(ml+yn*c+ab_margin/2-10+10),(ml+k*c+c),\
					(ml+yn*c+ab_margin/2+10+10)],fill=(rc,gc,bc),outline=(rc,gc,bc))
				k = k + 1

		draw_png.line([ml,(ml+yn*c+ab_margin/2+60+20),(ml+k*c),\
				(ml+yn*c+ab_margin/2+60+20)],fill=(64,64,64))
		draw_png.line([ml,(ml+yn*c+ab_margin/2-60),(ml+k*c),\
				(ml+yn*c+ab_margin/2-60)],fill=(64,64,64))
		draw_png.line([ml,(ml+yn*c+ab_margin/2-60),ml,\
				(ml+yn*c+ab_margin/2+60+20)],fill=(64,64,64))
		draw_png.line([(ml+k*c),(ml+yn*c+ab_margin/2-60),\
				(ml+k*c),(ml+yn*c+ab_margin/2+60+20)],fill=(64,64,64))

		draw_png.text([(ml+k*c+30),(ml+yn*c+ab_margin/2-20)],fill=(0,0,0),text=\
				("proportion of \"A\" alleles"),font=font)
		draw_png.text([(ml+k*c+30),(ml+yn*c+ab_margin/2+20)],fill=(0,0,0),text=\
				("proportion of \"B\" alleles"),font=font)

	#######################################
	if draw_ril_image == "TRUE":
		k = 0
		while k < xrn:
			# for ril_id in ril_list:
			z = k + 1
			ril_id = str(z)
			try:
				count_A  = rils_A[ril_id]
				count_B  = rils_B[ril_id]
				count_ALL = count_A + count_B
				# a_value  = int(round((count_A*1.0/count_ALL)*100))
				# b_value  = int(round((count_B*1.0/count_ALL)*100))
				ab_value_ril = int(round(((count_A-count_B)*1.0/count_ALL)*100))
				a_value  = int(round((count_A*0.75/count_ALL)*100))
				b_value  = int(round((count_B*0.75/count_ALL)*100))
				ab_value = int(round(((count_A-count_B)*0.75/count_ALL)*100))
				ab_diff  = a_value - b_value
				# print "COUNT: " + `count_A` + '\t' + `count_B` + '\t' + `count_ALL` + '\t' + ril_id
				# k = k + 1
			except:
				count_ALL = 0
				count_A = 0
				count_B = 0
				ab_value_ril = "-"
				# k = k + 1
			### RIL TEXT FILE
			try:
				value_r = double_ril[ril_id]
			except:
				value_r = 0
			ab_value_ril = str(ab_value_ril)
			out_file6.write(ril_id + '\t' + `count_A` + '\t' + `count_B` + '\t' + `count_ALL` \
				+ '\t' + ab_value_ril + '\t' + "---" + '\t' + `value_r` + '\n')
			### IMAGE FILE
			if count_ALL != 0:
				if ab_value >=0:
					rc = 255
					gc = 0
					bc = 0
				if ab_value < 0:
					rc = 0
					gc = 0
					bc = 255
				draw_ril.rectangle([(ml+k*c),(ml+yn*c+ab_margin/2+30),(ml+k*c+c),\
					(ml+yn*c+ab_margin/2-a_value+30)],fill=(200,0,0),outline=(200,0,0))
				draw_ril.rectangle([(ml+k*c),(ml+yn*c+ab_margin/2+30),(ml+k*c+c),\
					(ml+yn*c+ab_margin/2+b_value+30)],fill=(0,0,200),outline=(0,0,200))
				draw_ril.rectangle([(ml+k*c),(ml+yn*c+ab_margin/2+30),(ml+k*c+c),\
					(ml+yn*c+ab_margin/2+30-ab_value)],fill=(rc,gc,bc),outline=(rc,gc,bc))
				k = k + 1
			if count_ALL == 0:
				rc = 82
				gc = 82
				bc = 82
				draw_ril.rectangle([(ml+k*c),(ml+yn*c+ab_margin/2-10+30),(ml+k*c+c),\
					(ml+yn*c+ab_margin/2+10+30)],fill=(rc,gc,bc),outline=(rc,gc,bc))
				k = k + 1

		## CLOSE RIL TEXT FILE
		out_file6.close()
		## CLOSE MARKER TEXT FILE
		out_file5.close()

		draw_ril.line([ml,(ml+yn*c+ab_margin/2+60+40),(ml+k*c),\
				(ml+yn*c+ab_margin/2+60+40)],fill=(64,64,64))
		draw_ril.line([ml,(ml+yn*c+ab_margin/2-60+20),(ml+k*c),\
				(ml+yn*c+ab_margin/2-60+20)],fill=(64,64,64))
		draw_ril.line([ml,(ml+yn*c+ab_margin/2-60+20),ml,\
				(ml+yn*c+ab_margin/2+60+40)],fill=(64,64,64))
		draw_ril.line([(ml+k*c),(ml+yn*c+ab_margin/2-60+20),\
				(ml+k*c),(ml+yn*c+ab_margin/2+60+20+20)],fill=(64,64,64))

		draw_ril.text([(ml+k*c+30),(ml+yn*c+ab_margin/2-20+20)],fill=(0,0,0),text=\
				("proportion of \"A\" genotype"),font=font)
		draw_ril.text([(ml+k*c+30),(ml+yn*c+ab_margin/2+20+20)],fill=(0,0,0),text=\
				("proportion of \"B\" genotype"),font=font)

	##########################################

	draw_png.text([ml,(ml+yn*c+35+ab_margin)],fill=(0,0,0),text=("MAP: " + list_id),font=font)
	draw_png.text([ml+350,(ml+yn*c+35+ab_margin)],fill=(0,0,0),text=("MATRIX: " + in_name),font=font)
	link_limit = round(link_limit,2)
	link_limit = str(link_limit)
	draw_png.text([ml+750,(ml+yn*c+35+ab_margin)],fill=(0,0,0),text=("CUTOFF: " + link_limit),font=font)

	### DOUBLE RECOMBINATION FOR MARKERS ###
	if draw_ril_image == "TRUE":
		draw_ril.text([(ml+k*c+30),(ml+yn*c+ab_margin/2-60)],fill=(0,0,0),text=\
				("number of double cross overs"),font=font)
		draw_ril.text([ml,(ml+yn*c+40+ab_margin)],fill=(0,0,0),text=("MAP: " + list_id),font=font)
		draw_ril.text([ml+350,(ml+yn*c+40+ab_margin)],fill=(0,0,0),text=("LOCUS FILE: " + loc_file),font=font)

	### DOUBLE RECOMBINATION FOR RILs ###
		z = 1
		while z <= xrn:
			# for ril_id in ril_list:
			ril_id = str(z)
			tick_update = math.fmod(z,2)
			if tick_update == 0:
				level = 15
			if tick_update != 0:
				level = 0
			try:
				value = double_ril[ril_id]
			except:
				value = 0
			if value == 0:
				draw_ril.text([(ml+(z-1)*c),(ml+yn*c+4+level)],fill=(0,0,125),text=(`value`),font=font)
			if value == 1:
				draw_ril.text([(ml+(z-1)*c),(ml+yn*c+4+level)],fill=(0,60,180),text=(`value`),font=font)
			if value == 2:
				draw_ril.text([(ml+(z-1)*c),(ml+yn*c+4+level)],fill=(125,0,0),text=(`value`),font=font)
			if value > 2 and value < 10:
				draw_ril.text([(ml+(z-1)*c),(ml+yn*c+4+level)],fill=(200,0,0),text=(`value`),font=font)
			if value >= 10:
				draw_ril.text([(ml+(z-1)*c),(ml+yn*c+4+level)],fill=(225,0,0),text=("X"),font=font)
			z = z + 1

	in_file.close()
	id_file.close()
	out_file1.close()
	print ""
	print "PROCESSING LARGE SIZE"
	my_lovely_image.save(out_file2)

        ###########################################
        if draw_ril_image == "TRUE":
                print "PROCESSING RIL IMAGE"
                my_ril_image.save(out_ril_large_name, my_ril_image.format)
                mini_ril = Image.open(out_ril_large_name)
                mini_ril = mini_ril.resize((200,175),Image.ANTIALIAS)
                mini_ril = mini_ril.filter(ImageFilter.SHARPEN)
                mini_ril.save(out_ril_small_name, mini_ril.format)
        ###########################################

	#########   SCALING THE IMAGE   ###########
	im = Image.open(out_file2)
	width, height = im.size
	print "==========================================="
	print "IMAGE SIZE:  " + `width` + " X " + `height`
	print "==========================================="
	if cgpdb_style == "FALSE":
		out_medium_name = out_name + '.matrix2d.medium.png'
		out_5000_name = out_name + '.matrix2d.5000.png'
		out_2000_name = out_name + '.matrix2d.2000.png'
		out_small_name = out_name + '.matrix2d.small.png'
	### DUMMY CGPDB ##################################
	if cgpdb_style == "TRUE":
		out_medium_name = out_name + '.medium.png'
		out_5000_name = out_name + '.5000.png'
		out_2000_name = out_name + '.2000.png'
		out_small_name = out_name + '.xsmall.png'
	##################################################
	huge_image = "FALSE"
	if width > 1500 and height > 1200:
		print "PROCESSING MEDIUM SIZE"
		im_m = im.resize((1500,1200),Image.ANTIALIAS)
		im_m = im_m.filter(ImageFilter.SHARPEN)
		im_m.save(out_medium_name, im_m.format)
		huge_image = "TRUE"
	######################################################
	# if width > 2000 and height > 1500:
	#	print "PROCESSING 2000 PIXELS SIZE"
	#	im_2 = im.resize((2000,1500),Image.ANTIALIAS)
	#	im_2 = im_2.filter(ImageFilter.SHARPEN)
	#	im_2.save(out_2000_name, im_2.format)
	######################################################
	# if width >= 10000:
	#	print "PROCESSING 5000 PIXELS SIZE"
	#	im_5 = im.resize((5500,5000),Image.ANTIALIAS)
	#	im_5 = im_5.filter(ImageFilter.SHARPEN)
	#	im_5.save(out_5000_name, im_5.format)
	######################################################
	print "PROCESSING SMALL SIZE"
	if cgpdb_style == "FALSE":
		if huge_image == "TRUE":
			im_s = im_m.resize((200,175),Image.ANTIALIAS)
		if huge_image == "FALSE":
			im_s = im.resize((200,175),Image.ANTIALIAS)
	if cgpdb_style == "TRUE":
		if huge_image == "TRUE":
			im_s = im_m.resize((125,100),Image.ANTIALIAS)
		if huge_image == "FALSE":
			im_s = im.resize((125,100),Image.ANTIALIAS)
	print " HUGE IMAGE - " + huge_image
	im_s = im_s.filter(ImageFilter.SHARPEN)
	im_s.save(out_small_name, im_s.format)
	###########################################
	# if draw_ril_image == "TRUE":
	#	print "PROCESSING RIL IMAGE"
	#	my_ril_image.save(out_ril_large_name, my_ril_image.format)
	#	mini_ril = Image.open(out_ril_large_name)
	#	mini_ril = mini_ril.resize((200,175),Image.ANTIALIAS)
	#	mini_ril = mini_ril.filter(ImageFilter.SHARPEN)
	#	mini_ril.save(out_ril_small_name, mini_ril.format)
	###########################################
	print "===================================="
	print "             DONE!                  "
	print "===================================="

import sys
import math
import time
import Image
import ImageDraw
import ImageFont
import ImageFilter

if __name__ == "__main__":

	global cgpdb_style

	cgpdb_style = "FALSE"
	# cgpdb_style = "TRUE"

	if len(sys.argv) <= 11 or len(sys.argv) > 12:
		print ""
		print "Program usage: "
		print "[matrix_file] [map_file] [output_file] [frame_marker_list] [red_list] [loc_file] [REC/BIT/LOD] [GRAPH_OPT] [LINK_CUT] [LARGE/SMALL] [RIL_TYPE]"
		print "frame_marker_list is optional, if you do not have it just type X"
		print "red_list is a list of markers to highlight in red"
		print "red_list is optional, if you do not have it just type Y"
		print "loc_file is optional, if you do not have it just type Z"
		print "-------------------------------------------------------"
		print "GRAPH/NOGRAPH option - if 8-th option is \"GRAPH\" then CIRCULAR graph will be generated"
		print "LINK_CUT - linkage cutoff value on circular graph (0.9 - 0.8)"
		print "generation of circular graph needs lot of computer memory"
		print "do not use it for set of markers larger than 1000"
		print "-------------------------------------------------------"
		print "LARGE/SMALL option - how many pixels per cell. LARGE - 10 pixels, SMALL - 6 pixels"
		print "use SMALL option for large datasets - more than 1000 markers"
		print "-------------------------------------------------------"
		print "RIL type must be \"RIL\" or \"F2\""
		print ""
		sys.exit()
	if len(sys.argv) == 12:
		in_matrix    = sys.argv[1]
		in_map       = sys.argv[2]
		out_name     = sys.argv[3]
		purple_list  = sys.argv[4]
		red_list     = sys.argv[5]
		loc_file     = sys.argv[6]
		data_type    = sys.argv[7]
		graph_opt    = sys.argv[8]
		link_cut     = sys.argv[9]
		cell_size    = sys.argv[10]
		ril_type     = sys.argv[11]

		if data_type != "REC" and data_type != "BIT" and data_type != "LOD":
			print ""
			print "last value must be REC or BIT or LOD"
			print ""
			sys.exit()
		if data_type == "REC":
			column_n = 1
			link_limit = 0.4
			data_type = "REC"
			link_limit = float(link_limit)
		if data_type == "BIT":
			column_n = 2
			link_limit = 50
			data_type = "BIT"
			link_limit = int(link_limit)
		if data_type == "LOD":
			column_n = 2
			link_limit = 3
			data_type = "LOD"
			link_limit = int(link_limit)
		if cell_size != "LARGE" and cell_size != "SMALL":
			print ""
			print "last argument must be LARGE or SMALL"
			print ""
			sys.exit()
		if graph_opt != "GRAPH" and  graph_opt != "NOGRAPH":
			print ""
			print "circular graph argument must be GRAPH or NOGRAPH"
			print ""
			sys.exit()
		if ril_type != "RIL" and  ril_type != "F2":
			print ""
			print "RIL type must be \"RIL\" or \"F2\""
			print ""
			sys.exit()
		
		diag_fill = 0
		dummy_fill = "-"
		# rgb_coeff and chick_sat in use only for REC scores
		rgb_coeff = 255
		chick_sat = 4
		column_n   = int(column_n)
		rgb_coeff  = float(rgb_coeff)
		chick_sat  = int(chick_sat)
		max_ril_n    = 0
		max_ril_n  = int(max_ril_n)
		# link_limit = float(link_limit)
		link_cut = float(link_cut)

		Seqs_Matrix(in_matrix, in_map, out_name, column_n, diag_fill, dummy_fill, \
			rgb_coeff, chick_sat, link_limit, purple_list, red_list, loc_file, data_type, max_ril_n, link_cut, cell_size, graph_opt, ril_type)
### THE END ###


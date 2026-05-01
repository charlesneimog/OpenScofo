; ================================
; COMMENTS
; ================================
(comment) @comment

((comment) @comment.documentation
 (#match? @comment "^/\\*\\*"))

; ================================
; STRUCTURE
; ================================
(CONFIG) @keyword.directive
(EVENT) @keyword

; ================================
; CONFIG
; ================================
(CONFIG
  value: (_) @number)

; ================================
; EVENTS (semantic types)
; ================================
(note_event) @type.builtin
(rest_event) @type.builtin
(ptech_event) @type.builtin
(utech_event) @type.builtin
(chord_event) @type.builtin
(trill_event) @type.builtin
(lua_event) @type.builtin




; ================================
; EVENT FIELDS
; ================================
(pitch) @string
[
  (utech_event
    technique: (identifier) @string)
  (ptech_event
    technique: (identifier) @string)
]

[
  (note_event
    duration: (number) @variable.parameter)
  (rest_event
    duration: (number) @variable.parameter)
  (ptech_event
    duration: (number) @variable.parameter)
  (utech_event
    duration: (number) @variable.parameter)
]

; ================================
; ACTIONS (runtime layer)
; ================================
(action) @function.call

(exec
  "sendto" @function)

(exec
  "luacall" @function)

(delay) @keyword.operator
(time_unit) @type

; receiver (e.g. pitchshift, delay, freeze)
(exec
  receiver: (identifier) @variable)

; ================================
; PD ARGUMENTS (structured)
; ================================
(pdargs) @punctuation.bracket

(pdarg) @number

(pdarg
  (identifier) @variable.parameter)

; ================================
; LUA EMBED
; ================================
(LUA) @keyword
(lua_body) @none
(lua_call) @function.call
(lua_comment) @comment

; ================================
; MISC
; ================================
(time_unit) @type
